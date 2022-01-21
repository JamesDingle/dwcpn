use crate::dwcpn::modules::config::TIMESTEPS;
use std::f64::consts::PI;

pub fn generate_zenith_array(
    time_array: [f64; TIMESTEPS],
    delta: f64,
    phi: f64,
) -> ([f64; TIMESTEPS], [f64; TIMESTEPS]) {
    let mut zenith_array: [f64; TIMESTEPS] = [0.0; TIMESTEPS];
    let mut zenith_d_array: [f64; TIMESTEPS] = [0.0; TIMESTEPS];

    for i in 0..TIMESTEPS {
        zenith_array[i] = compute_zenith(time_array[i], delta, phi);
        zenith_d_array[i] = zenith_array[i] * (180.0 / PI);
    }

    (zenith_array, zenith_d_array)
}

pub fn compute_zenith(local_time: f64, delta: f64, phi: f64) -> f64 {
    let th: f64 = (local_time - 12.0) * (PI / 12.);

    let mut zen: f64 = delta.sin() * phi.sin() + delta.cos() * phi.cos() * th.cos();

    if zen < -1.0 {
        zen = -1.0;
    } else if zen > 1.0 {
        zen = 1.0;
    }

    zen = (PI / 2.0) - zen.asin();

    return zen;
}

pub fn compute_zenith_time(delta: f64, phi: f64, zenith_d: f64) -> f64 {
    let zenith_r = zenith_d * PI / 180.0;

    let cos_theta = (zenith_r.cos() - (delta.sin() * phi.sin())) / (delta.cos() * phi.cos());

    let theta_1 = cos_theta.acos();

    let mut local_time: f64;

    if cos_theta <= -1.0 {
        local_time = 0.0;
    } else if cos_theta >= 1.0 {
        local_time = -1.0;
    } else {
        local_time = 12.0 + (theta_1 / (PI / 12.0));
    }

    if local_time > 12.0 {
        local_time = 24.0 - local_time;
    }

    local_time
}