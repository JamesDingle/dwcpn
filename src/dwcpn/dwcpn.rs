use crate::dwcpn::modules::chl_profile::{gen_chl_profile};
use crate::dwcpn::modules::config::{DEPTH_PROFILE_STEP, TIMESTEPS};
use crate::dwcpn::modules::irradiance::{compute_irradiance_components, correct_and_recompute_irradiance_components, lookup_thekaekara_correction};
use crate::dwcpn::modules::pp_profile::{compute_pp_depth_profile, compute_prochloro_profile};
use crate::dwcpn::modules::time::{compute_sunrise, generate_time_array};
use crate::dwcpn::modules::zenith::{generate_zenith_array, compute_zenith_time};
use std::f64::consts::PI;
use crate::dwcpn::modules::light_profile::calc_light_decay_profile;
use crate::{DEPTH_PROFILE_COUNT, ModelInputs, ModelOutputs, ModelSettings, PPErrors};


pub fn calc_production(input: &ModelInputs, settings: &ModelSettings) -> Result<ModelOutputs, PPErrors> {

    // generate chl depth profile
    let (depth_array, chl_profile) = gen_chl_profile(input, settings);

    // compute sunrise and generate time array
    let (sunrise, delta, phi) = compute_sunrise(input.iday, input.lat);

    let zenith_80_time = compute_zenith_time(delta, phi, 80.0);

    let (time_array, delta_t) = generate_time_array(zenith_80_time);
    let (zenith_array, zenith_d_array) = generate_zenith_array(time_array, delta, phi);

    let mut start_time_idx: f64 = -1.0;
    let mut day_length: f64 = 0.0;
    let mut iom: f64 = 0.0;
    let mut delta_prestart: f64 = 0.0;

    let solar_correction = lookup_thekaekara_correction(input.iday);

    // arrays to store results
    let mut pp: [f64; TIMESTEPS] = [0.0; TIMESTEPS];
    let mut euphotic_depth: [f64; TIMESTEPS] = [0.0; TIMESTEPS];

    if settings.prochloro_inputs.is_some() {
        let mut pro_1_profile: Option<[f64; DEPTH_PROFILE_COUNT]> = Some([0.0; DEPTH_PROFILE_COUNT]);
        let mut pro_2_profile: Option<[f64; DEPTH_PROFILE_COUNT]> = Some([0.0; DEPTH_PROFILE_COUNT]);
        let mut pro_total_profile: Option<[f64; DEPTH_PROFILE_COUNT]> = Some([0.0; DEPTH_PROFILE_COUNT]);
    } else {
        let mut pro_1_profile: Option<[f64; DEPTH_PROFILE_COUNT]> = None;
        let mut pro_2_profile: Option<[f64; DEPTH_PROFILE_COUNT]> = None;
        let mut pro_total_profile: Option<[f64; DEPTH_PROFILE_COUNT]> = None;
    }
    let mut pro_total_count: usize = 0;

    // spectral i star is calculated as a running mean
    let mut spectral_i_star_sum: f64 = 0.0;
    let mut spectral_i_star_count: f64 = 0.0;

    let mut start_time: f64;

    // loop over time array (from sunrise to noon)
    for t in 0..TIMESTEPS {
        // if the zenith angle is yet to go below 80° then skip to the next time step
        if zenith_d_array[t] >= 80.00005 {
            continue;
        }

        // update start_time if necessary and assign delta_prestart variable
        // to be used for daily integration purposes later.
        // t_start is t_idx value when Zenith angle becomes >80 degrees
        // start_time is calculation start time, start_time_index is the index
        // delta_prestart is time elapsed between dawn and start_time
        if start_time_idx < 0.0 {
            start_time_idx = t as f64;

            // this line has start_time_idx corrected with +1 as the original code
            // came from FORTRAN which uses indices that start at 1
            // start_time = sunrise.clone() + delta_t.clone() * (start_time_idx);
            start_time = zenith_80_time;

            day_length = 2.0 * (12.0 - sunrise);

            // iom = noon time maximum
            iom = input.par * PI / (2.0 * day_length);
            delta_prestart = start_time - sunrise;

            if settings.iom_only == true {
                return Ok(
                    ModelOutputs {
                        pp_day: None,
                        euphotic_depth: None,
                        spectral_i_star: None,
                        par_noon_max: Some(iom),
                        pro_1_profile: None,
                        pro_2_profile: None,
                        pro_total_profile: None
                    }
                )
            }

        }

        // compute direct and diffuse irradiance components at sea level
        let (direct, diffuse) =
            compute_irradiance_components(zenith_array[t], zenith_d_array[t]);

        let (direct_corrected, diffuse_corrected) = correct_and_recompute_irradiance_components(
            direct,
            diffuse,
            solar_correction,
            iom,
            day_length,
            sunrise,
            zenith_array[t],
            time_array[t],
            input.cloud
        );

        let (i_alpha_profile, par_profile) = calc_light_decay_profile(
            chl_profile,
            direct_corrected,
            diffuse_corrected,
            zenith_array[t],
            &input
        );


        let pp_profile = compute_pp_depth_profile(
            &chl_profile,
            &depth_array,
            &i_alpha_profile,
            &par_profile,
            &input
        );

        match pp_profile {
            Ok(mut pp_profile) => {
                euphotic_depth[t] = pp_profile.euphotic_depth;

                if pp_profile.euph_index == 0 { pp_profile.euph_index = 1; }

                for z in 0..pp_profile.euph_index {
                    pp[t] = pp[t] + DEPTH_PROFILE_STEP * (pp_profile.pp_profile[z] + pp_profile.pp_profile[z + 1]) / 2.0;
                }

                pp[t] = pp[t]
                    + pp_profile.pp_profile[pp_profile.euph_index]
                    * (euphotic_depth[t]
                    - (pp_profile.euph_index as f64 - 1.0) * DEPTH_PROFILE_STEP);

                spectral_i_star_sum = spectral_i_star_sum + (pp_profile.spectral_i_star / (pp_profile.euph_index as f64).abs());
                spectral_i_star_count = spectral_i_star_count + 1.0;
            },
            Err(e) => println!("{:?}", e)
        }


        if settings.prochloro_inputs.is_some() {
            let pro_surf = &settings.prochloro_inputs.as_ref().unwrap().prochloro_surface;
            let pro_max = &settings.prochloro_inputs.as_ref().unwrap().prochloro_maximum;

            let prochloro_profile = compute_prochloro_profile(
                &chl_profile,
                &depth_array,
                &i_alpha_profile,
                &par_profile,
                &input,
                &pro_surf,
                &pro_max
            );

            match prochloro_profile {
                Ok(mut prochloro_profile) => {
                    euphotic_depth[t] = prochloro_profile.euphotic_depth;

                    if prochloro_profile.euph_index == 0 { prochloro_profile.euph_index = 1; }

                    let mut pro_output_profile: [f64; DEPTH_PROFILE_COUNT] = pro_total_profile.unwrap();
                    let mut pro_1_output_profile: [f64; DEPTH_PROFILE_COUNT] = pro_1_profile.unwrap();
                    let mut pro_2_output_profile: [f64; DEPTH_PROFILE_COUNT] = pro_2_profile.unwrap();

                    // TODO: double check if we need to cut off at euphotic depth for prochlorococcus
                    for z in 0..prochloro_profile.euph_index {
                        pro_output_profile[z] = pro_output_profile[z] + prochloro_profile.pro_1_profile[z] + prochloro_profile.pro_2_profile[z];
                        pro_1_output_profile[z] = pro_1_output_profile[z] + prochloro_profile.pro_1_profile[z];
                        pro_2_output_profile[z] = pro_2_output_profile[z] + prochloro_profile.pro_2_profile[z];
                    }

                    pro_total_count += 1;
                    pro_total_profile = Some(pro_output_profile);
                    pro_1_profile = Some(pro_1_output_profile);
                    pro_2_profile = Some(pro_2_output_profile);
                },
                Err(e) => println!("{:?}", e)
            }
        }

    } // time loop

    let mut pp_day = pp[0] * delta_prestart / 2.0;
    let mut max_euphotic_depth: f64 = 0.0;
    for t in 0..TIMESTEPS - 1 {
        pp_day = pp_day + ((pp[t] + pp[t + 1]) * delta_t / 2.0);

        if max_euphotic_depth.abs() < euphotic_depth[t].abs() {
            max_euphotic_depth = euphotic_depth[t].abs();
        }
    }


    // // calculate final mean for spectral i star
    let mut spectral_i_star_mean: f64 = 0.0;
    if spectral_i_star_count > 0.0 {
        spectral_i_star_mean = spectral_i_star_sum / spectral_i_star_count;
    }

    // mutliply by two because we have only integrated over half of the day
    pp_day = pp_day * 2.0;

    // Calculate mean (along time) prochlorococcus for every depth
    if pro_total_profile.is_some() {
        for z in 0..DEPTH_PROFILE_COUNT {
            pro_total_profile.unwrap()[z] /= pro_total_count as f64;
            pro_1_profile.unwrap()[z] /= pro_total_count as f64;
            pro_2_profile.unwrap()[z] /= pro_total_count as f64;
        }
    }

    if pp_day > 10000.0 {
        Err(PPErrors::PPTooHigh)
    } else {
        Ok(
            ModelOutputs {
                pp_day: Some(pp_day),
                euphotic_depth: Some(max_euphotic_depth),
                spectral_i_star: Some(spectral_i_star_mean),
                par_noon_max: Some(iom),
                pro_1_profile,
                pro_2_profile,
                pro_total_profile
            }
        )
    }


}
