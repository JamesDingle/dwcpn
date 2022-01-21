use crate::dwcpn::modules::absorption::calc_ac;
use crate::dwcpn::modules::config::{AW, DEPTH_PROFILE_COUNT, DEPTH_PROFILE_STEP, WL_ARRAY, WL_COUNT, DELTA_LAMBDA};
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

    let mut i_zero: [f64; WL_COUNT] = [0.0; WL_COUNT];
    let mut mu_d: [f64; WL_COUNT] = [0.0; WL_COUNT];
    let mut i_z: [f64; WL_COUNT] = [0.0; WL_COUNT];
    let mut k: [f64; WL_COUNT] = [0.0; WL_COUNT];

    let zenith_w: f64 = (zenith_r.sin() / 1.333).asin();

    for l in 0..WL_COUNT {
        i_zero[l] = direct_irradiance[l] + diffuse_irradiance[l];
        mu_d[l] =
            (direct_irradiance[l] * zenith_w.cos() + diffuse_irradiance[l] * 0.831000) / i_zero[l];
        i_z[l] = i_zero[l];
    }

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

        // let ac440 = ac[8]; // TODO: fix this to search for/interpolate to 440nm
        let ac440 = linear_interp(&WL_ARRAY, &ac, 440.0);

        let power = -(chl.log10());
        let ay440 = yellow_substance * ac440;

        let bc660 = 0.407 * chl.powf(0.795);
        let mut bbtilda = (0.78 + 0.42 * power) * 0.01;

        if bbtilda < 0.0005 {
            bbtilda = 0.0005;
        } else if bbtilda > 0.01 {
            bbtilda = 0.01;
        }

        for l in 0..WL_COUNT {
            let wl = WL_ARRAY[l];
            let a = AW[l] + ac[l] + ay440 * ay[l] + 2.0 * bbr[l];
            let mut bc = bc660 * (660.0 / wl).powf(power);

            if bc < 0.0 {
                bc = 0.0
            };

            let bb = bc * bbtilda + bw[l] * 0.50;

            k[l] = (a + bb) / mu_d[l];

            par_profile[z] = par_profile[z] + i_z[l] * DELTA_LAMBDA;
        }

        let mut i_alpha = 0.0;
        for l in 0..WL_COUNT {

            // this conversion expects pi_alpha to be in units of
            // mgC mgChl^-1 h^-1 (W m^-2)^-1
            // a.ka. (mgC per mgChl per Hour) / (Watts per m^2)
            // the line below converts irradiance (light units) to einsteins per m^2 per hour
            // this makes it compatible with the par units
            let x = province_alpha * ac[l] * 6022.0 / (2.77 * 36.0 * ac_mean);


            i_alpha = i_alpha + x * DELTA_LAMBDA * i_z[l] / mu_d[l];
            i_z[l] = i_z[l] * (-k[l] * DEPTH_PROFILE_STEP).exp();
        }


        // pp equation has been updated after discussion with Shubha 2018/08/16
        // pp_profile[z] = (i_alpha / (1.0 + (i_alpha / province_pmb).powf(2.0)).sqrt()) * chl; // this is the old primary production equation.

        pp_profile[z] = chl * province_pmb * (1.0 - (-i_alpha / province_pmb).exp());
        // spectral_i_star_profile[z] = i_alpha / province_pmb.clone();
        i_alpha_sum = i_alpha_sum + i_alpha;


        if z > 0 {
            if par_profile[z] < (0.01 * par_profile[0]) {
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
