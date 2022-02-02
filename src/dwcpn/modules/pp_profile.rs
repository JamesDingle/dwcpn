use crate::dwcpn::modules::absorption::calc_ac;
use crate::dwcpn::modules::config::{AW, DEPTH_PROFILE_COUNT, DEPTH_PROFILE_STEP, WL_ARRAY, WL_COUNT, DELTA_LAMBDA};
use crate::dwcpn::modules::light_profile::{calc_i_z_decay, calc_par_profile, init_mu_d_and_i_z};
use crate::dwcpn::modules::linear_interp::linear_interp;


pub struct PpProfile {
    pub pp_profile: [f64; DEPTH_PROFILE_COUNT],
    pub par_profile: [f64; DEPTH_PROFILE_COUNT],
    pub euphotic_depth: f64,
    pub euph_index: usize,
    pub spectral_i_star: f64,
    pub success: bool,
}

pub fn calculate_bw() -> [f64; WL_COUNT] {
    // scattering coefficient of pure seawater at 500nm
    const BW500: f64 = 0.00288;

    let mut bw: [f64; WL_COUNT] = [0.0; WL_COUNT];

    for i in 0..WL_COUNT {
        bw[i] = BW500 * (WL_ARRAY[i] / 500.0).powf(-4.3);
    }

    bw
}

pub fn calculate_bbr() -> [f64; WL_COUNT] {
    const BR488: f64 = 0.00027;

    let mut bbr: [f64; WL_COUNT] = [0.0; WL_COUNT];

    for i in 0..WL_COUNT {
        bbr[i] = 0.5 * BR488 * (WL_ARRAY[i] / 488.0).powf(-5.3);
    }

    bbr
}

pub fn calculate_ay() -> [f64; WL_COUNT] {
    let mut ay: [f64; WL_COUNT] = [0.0; WL_COUNT];

    for i in 0..WL_COUNT {
        ay[i] = (-0.014 * (WL_ARRAY[i] - 440.0)).exp();
    }

    ay
}

pub fn compute_pp_depth_profile(
    chl_profile: [f64; DEPTH_PROFILE_COUNT],
    depth_profile: [f64; DEPTH_PROFILE_COUNT],
    zenith_r: f64,
    direct_irradiance: [f64; WL_COUNT],
    diffuse_irradiance: [f64; WL_COUNT],
    bw: [f64; WL_COUNT],
    bbr: [f64; WL_COUNT],
    ay: [f64; WL_COUNT],
    province_alpha: f64,
    province_pmb: f64,
    yellow_substance: f64,
) -> PpProfile {
    let mut pp_profile: [f64; DEPTH_PROFILE_COUNT] = [0.0; DEPTH_PROFILE_COUNT];
    let mut par_profile: [f64; DEPTH_PROFILE_COUNT] = [0.0; DEPTH_PROFILE_COUNT];
    let mut euphotic_depth: f64 = 0.0;
    let mut success = false;

    let (mu_d, mut i_z) = init_mu_d_and_i_z(direct_irradiance, diffuse_irradiance, zenith_r);

    let mut i_alpha_sum: f64 = 0.0;
    let mut spectral_i_star: f64 = 0.0;

    for z in 0..DEPTH_PROFILE_COUNT {
        let chl = chl_profile[z];
        let (ac, ac_mean) = calc_ac(chl);

        if ac_mean == 0.0 {
            if z > 1 {
                let euph_index = z - 1;
                euphotic_depth = depth_profile[euph_index]
                    + DEPTH_PROFILE_STEP * (100.0 * par_profile[euph_index] / par_profile[0]).ln()
                    / (par_profile[euph_index] / par_profile[z]).ln();
                success = true;
                spectral_i_star = i_alpha_sum / province_pmb;
                return PpProfile {
                    pp_profile,
                    par_profile,
                    euphotic_depth,
                    euph_index,
                    spectral_i_star,
                    success,
                };
            }
        }

        //use of _temp variables to assign to outer variables (i_z and par_profile) are due to not
        //knowing how to correctly re-assign multiple values from a function
        // i_z is a running value of remaining light at each wavelength
        // par_z is the PAR available at current depth (z) - calculated from i_z
        let (i_alpha_test, i_z_temp, par_z) = calc_i_z_decay(
            ac,
            mu_d,
            i_z,
            chl,
            yellow_substance,
            ay,
            bbr,
            bw,
            province_alpha,
            ac_mean
        );

        i_z = i_z_temp; // stored to pass back into the function
        par_profile[z] = par_z;

        {
            // pp equation has been updated after discussion with Shubha 2018/08/16
            // pp_profile[z] = (i_alpha / (1.0 + (i_alpha / province_pmb).powf(2.0)).sqrt()) * chl; // this is the old primary production equation.

            // spectral_i_star_profile[z] = i_alpha / province_pmb.clone();
        }

        pp_profile[z] = chl * province_pmb * (1.0 - (-i_alpha_test / province_pmb).exp());
        i_alpha_sum = i_alpha_sum + i_alpha_test;

        if z > 0 && par_profile[z] < (0.01 * par_profile[0]) {
            let (euph_index, euphotic_depth) = integrate_euphotic_depth(z, depth_profile, par_profile);
            return PpProfile {
                pp_profile,
                par_profile,
                euphotic_depth,
                euph_index,
                spectral_i_star: i_alpha_sum / province_pmb,
                success: true,
            };
        }
    } // depth loop

    let euph_index = 0;
    return PpProfile {
        pp_profile,
        par_profile,
        euphotic_depth,
        euph_index,
        spectral_i_star,
        success,
    };
}

fn integrate_euphotic_depth(
    depth_index: usize,
    depth_profile: [f64; DEPTH_PROFILE_COUNT],
    par_profile: [f64; DEPTH_PROFILE_COUNT]
) -> (usize, f64) {
    let euph_index = depth_index - 1;
    let euphotic_depth = depth_profile[euph_index]
        + DEPTH_PROFILE_STEP * (100.0 * par_profile[euph_index] / par_profile[0]).ln()
        / (par_profile[euph_index] / par_profile[depth_index]).ln();

    (euph_index, euphotic_depth)
}