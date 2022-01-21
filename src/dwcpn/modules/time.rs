use crate::dwcpn::modules::config::TIMESTEPS;

pub fn generate_time_array(sunrise: f64) -> ([f64; TIMESTEPS], f64) {
    let mut time_array: [f64; TIMESTEPS] = [0.0; TIMESTEPS];
    let count = (TIMESTEPS) as f64;
    let delta_t: f64 = (12.0 - sunrise) / (count - 1.0);

    for i in 0..TIMESTEPS {
        let i_f64 = i as f64;
        time_array[i] = sunrise + delta_t * i_f64;
    }

    return (time_array, delta_t);
}

// pub fn generate_time_array_zenith(zenith_time: f64) -> ([f64; TIMESTEPS], f64) {
//     let mut time_array: [f64; TIMESTEPS] = [0.0; TIMESTEPS];
//     let count = (TIMESTEPS + 1) as f64;
//     let delta_t: f64 = (12.0 - zenith_time) / count;
//
//     for i in 0..TIMESTEPS {
//         let i_f64 = i as f64;
//         time_array[i] = zenith_time + delta_t * i_f64;
//     }
//
//     return (time_array, delta_t);
// }



pub fn compute_sunrise(jday: u16, lat: f64) -> (f64, f64, f64) {
    // convert julian day to double for the equation
    let dayfloat = jday as f64;

    let tau: f64 = std::f64::consts::PI * 2.00f64;

    let theta: f64 = tau * dayfloat / 365.0f64;

    let delta: f64 = 0.006918 - 0.399912 * theta.cos()
        + 0.070257 * theta.sin()
        - 0.006758 * (2.0 * theta).cos()
        + 0.000907 * (2.0 * theta).sin()
        - 0.002697 * (3.0 * theta).cos()
        + 0.001480 * (3.0 * theta).sin();

    let phi: f64 = lat * (tau / 360.0);

    let mut phidel: f64 = -phi.tan() * delta.tan();

    if phidel < -1.0 {
        phidel = -1.0; // 24hr sunlight
    } else if phidel > 1.0 {
        phidel = 1.0; // 24hr darkness
    }

    // return sunrise
    let sunrise: f64 = 12.0 - phidel.acos() * (360.0 / tau) / 15.0;

    return (sunrise, delta, phi);
}
