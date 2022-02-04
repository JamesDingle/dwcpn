use crate::dwcpn::modules::absorption::calc_ac;
use crate::dwcpn::modules::config::{AW, DELTA_LAMBDA, DEPTH_PROFILE_COUNT, DEPTH_PROFILE_STEP, WL_ARRAY, WL_COUNT};
use crate::dwcpn::modules::linear_interp::linear_interp;

pub struct LightProfile {
    pub par_profile: [f64; DEPTH_PROFILE_COUNT],
    pub i_alpha_profile: [f64; DEPTH_PROFILE_COUNT]
}

pub fn init_mu_d_and_i_z (
    direct_irradiance: [f64; WL_COUNT],
    diffuse_irradiance: [f64; WL_COUNT],
    zenith_r: f64
) -> ([f64; WL_COUNT], [f64; WL_COUNT]) {
    let mut i_zero: [f64; WL_COUNT] = [0.0; WL_COUNT];
    let mut mu_d: [f64; WL_COUNT] = [0.0; WL_COUNT];
    let mut i_z: [f64; WL_COUNT] = [0.0; WL_COUNT];
    let zenith_w: f64 = (zenith_r.sin() / 1.333).asin();

    for l in 0..WL_COUNT {
        i_zero[l] = direct_irradiance[l] + diffuse_irradiance[l];
        mu_d[l] =
            (direct_irradiance[l] * zenith_w.cos() + diffuse_irradiance[l] * 0.831000) / i_zero[l];
        i_z[l] = i_zero[l];
    }

    return (mu_d, i_z)
}

pub fn calc_i_z_decay(
    ac: [f64; WL_COUNT],
    mu_d: [f64; WL_COUNT],
    i_z: [f64; WL_COUNT],
    chl: f64,
    yellow_substance: f64,
    ay: [f64; WL_COUNT],
    bbr: [f64; WL_COUNT],
    bw: [f64; WL_COUNT],
    province_alpha: f64,
    ac_mean: f64
) -> (f64, [f64; WL_COUNT], f64) {
    let mut i_z= i_z.clone();
    let mut i_alpha = 0.0;
    let mut k: [f64; WL_COUNT] = [0.0; WL_COUNT];

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

    let mut par = 0.0;

    for l in 0..WL_COUNT {
        par = par + i_z[l] * DELTA_LAMBDA;
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

        // par_profile[z] = par_profile[z] + i_z[l] * DELTA_LAMBDA;

        // this conversion expects pi_alpha to be in units of
        // mgC mgChl^-1 h^-1 (W m^-2)^-1
        // a.ka. (mgC per mgChl per Hour) / (Watts per m^2)
        // the line below converts irradiance (light units) to einsteins per m^2 per hour
        // this makes it compatible with the par units
        let x = province_alpha * ac[l] * 6022.0 / (2.77 * 36.0 * ac_mean);



        i_alpha = i_alpha + x * DELTA_LAMBDA * i_z[l] / mu_d[l];
        i_z[l] = i_z[l] * (-k[l] * DEPTH_PROFILE_STEP).exp();
    }

// note to self: implement the correct passing of par_profile
    (i_alpha, i_z, par)
}

pub fn calc_light_decay_profile(
    chl_profile: [f64; DEPTH_PROFILE_COUNT],
    direct_irradiance: [f64; WL_COUNT],
    diffuse_irradiance: [f64; WL_COUNT],
    zenith_r: f64,
    yellow_substance: f64,
    ay: [f64; WL_COUNT],
    bbr: [f64; WL_COUNT],
    bw: [f64; WL_COUNT],
    province_alpha: f64
) -> ([f64; DEPTH_PROFILE_COUNT], [f64; DEPTH_PROFILE_COUNT]) {
    let mut i_alpha_profile = [0.0; DEPTH_PROFILE_COUNT];
    let mut par_profile = [0.0; DEPTH_PROFILE_COUNT];

    let (mu_d, mut i_z) = init_mu_d_and_i_z(direct_irradiance, diffuse_irradiance, zenith_r);

    for z in 0..DEPTH_PROFILE_COUNT {
        let chl = chl_profile[z];
        let (ac, ac_mean) = calc_ac(chl);

        if ac_mean == 0.0 { break; }

        let (i_alpha_z, i_z_temp, par_z) = calc_i_z_decay(
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
        i_alpha_profile[z] = i_alpha_z;
        i_z = i_z_temp;
        par_profile[z] = par_z;
    }

    (i_alpha_profile,par_profile)
}

pub fn calc_par_profile(i_z: [f64; WL_COUNT]) -> [f64; DEPTH_PROFILE_COUNT] {
    let mut par_profile = [0.0; DEPTH_PROFILE_COUNT];

    for z in 0..DEPTH_PROFILE_COUNT {
        for l in 0..WL_COUNT {
            par_profile[z] = par_profile[z] + i_z[l] * DELTA_LAMBDA;
        }
    }

    par_profile
}