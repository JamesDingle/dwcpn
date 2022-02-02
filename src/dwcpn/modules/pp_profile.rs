use crate::dwcpn::dwcpn::{ModelInputs, PPErrors};
use crate::dwcpn::dwcpn::PPErrors::DWCPNError;
use crate::dwcpn::modules::absorption::calc_ac;
use crate::dwcpn::modules::config::{AW, DEPTH_PROFILE_COUNT, DEPTH_PROFILE_STEP, WL_ARRAY, WL_COUNT, DELTA_LAMBDA};
use crate::dwcpn::modules::light_profile::{calc_i_z_decay, calc_light_decay_profile, calc_par_profile, init_mu_d_and_i_z};
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
    i_alpha_profile: [f64; DEPTH_PROFILE_COUNT],
    par_profile: [f64; DEPTH_PROFILE_COUNT],
    model_inputs: &ModelInputs
) -> Result<PpProfile, PPErrors> {
    let mut pp_profile: [f64; DEPTH_PROFILE_COUNT] = [0.0; DEPTH_PROFILE_COUNT];

    let mut i_alpha_sum: f64 = 0.0;

    for z in 0..DEPTH_PROFILE_COUNT {
        pp_profile[z] = chl_profile[z] * model_inputs.pmb * (1.0 - (-i_alpha_profile[z] / model_inputs.pmb).exp());
        i_alpha_sum = i_alpha_sum + i_alpha_profile[z];

        if z > 0 && par_profile[z] < (0.01 * par_profile[0]) {
            let (mut euph_index,  mut euphotic_depth) = integrate_euphotic_depth(z, depth_profile, par_profile);

            // clamp euphotic_depth to physical depth of ocean if it is lower
            if euphotic_depth.abs() > model_inputs.z_bottom.abs() {
                euph_index = DEPTH_PROFILE_COUNT - 1;
                euphotic_depth = model_inputs.z_bottom.abs();
            }

            return Ok(PpProfile {
                pp_profile,
                par_profile,
                euphotic_depth,
                euph_index,
                spectral_i_star: i_alpha_sum / model_inputs.pmb,
                success: true,
            });
        }
    } // depth loop

    return Err(DWCPNError);
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