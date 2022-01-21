use crate::dwcpn::modules::config::{WL_ARRAY, WL_COUNT};
use crate::dwcpn::modules::linear_interp::linear_interp;

// DO NOT CHANGE THIS UNLESS YOU HAVE NEW LOOKUP TABLES FOR ALL OF THE BELOW CONST ARRAYS
const TRANSMITTANCE_WL_COUNT: usize = 24;
const TRANSMITTANCE_WAVELENGTHS: [f64; TRANSMITTANCE_WL_COUNT] = [
    400.000, 410.000, 420.000, 430.000, 440.000, 450.000, 460.000, 470.000, 480.000, 490.000,
    500.000, 510.000, 520.000, 530.000, 540.000, 550.000, 570.000, 593.000, 610.000, 630.000,
    656.000, 667.600, 690.000, 710.000,
];

// ozone absorption coefficients
const OZONE_ABS: [f64; TRANSMITTANCE_WL_COUNT] = [
    0.000, 0.000, 0.000, 0.000, 0.000, 0.003, 0.006, 0.009, 0.014, 0.021, 0.030, 0.040, 0.048,
    0.063, 0.075, 0.095, 0.120, 0.119, 0.132, 0.120, 0.065, 0.060, 0.028, 0.018,
];

// Water vapour absorption is set to 2.0
const W: f64 = 2.0;
// water vapour absorption coefficients
const WATER_VAPOUR_ABS: [f64; TRANSMITTANCE_WL_COUNT] = [
    0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
    0.000, 0.000, 0.000, 0.000, 0.075, 0.000, 0.000, 0.000, 0.000, 0.016, 0.0125,
];

// Aerosol transmittance coefficients
const ALPHA1: f64 = 1.0274;
const BETA1: f64 = 0.1324;
const ALPHA2: f64 = 1.206;
const BETA2: f64 = 0.117;

// correct factors for diffuse irradiance for 5 wavelengths at 7 zenith angles
const CORRECTION: [[f64; 5]; 7] = [
    [1.11, 1.04, 1.15, 1.12, 1.32],
    [1.13, 1.05, 1.00, 0.96, 1.12],
    [1.18, 1.09, 1.00, 0.96, 1.07],
    [1.24, 1.11, 0.99, 0.94, 1.02],
    [1.46, 1.24, 1.06, 0.99, 1.10],
    [1.70, 1.34, 1.07, 0.96, 0.90],
    [2.61, 1.72, 1.22, 1.04, 0.80],
];

const CORRECTION_ZEN_LOOKUP: [f64; 7] = [0., 37., 48.19, 60., 70., 75., 80.];

// extra-terrestrial spectral irradiance
const ET_SPECTRAL_IRRADIANCE: [f64; TRANSMITTANCE_WL_COUNT] = [
    1479.1, 1701.3, 1740.4, 1587.2, 1837.0, 2005.0, 2043.0, 1987.0, 2027.0, 1896.0, 1909.0, 1927.0,
    1831.0, 1891.0, 1898.0, 1892.0, 1840.0, 1768.0, 1728.0, 1658.0, 1524.0, 1531.0, 1420.0, 1399.0,
];

pub fn lookup_thekaekara_correction(julian_day: u16) -> f64 {
    let day_points: [f64; 25] = [
        0., 3., 31., 42., 59., 78., 90., 93., 120., 133., 151., 170., 181., 183., 206., 212., 243.,
        265., 273., 277., 304., 306., 334., 355., 365.,
    ];

    let irradiance_points: [f64; 25] = [
        1399., 1399., 1393., 1389., 1378., 1364., 1355., 1353., 1332., 1324., 1316., 1310., 1309.,
        1309., 1312., 1313., 1329., 1344., 1350., 1353., 1347., 1375., 1392., 1398., 1399.,
    ];

    let mut idx: usize = 0;
    for i in 0..25 {
        idx = i;
        if i >= julian_day as usize {
            break;
        }
    }

    if idx == 0 {
        return irradiance_points[0];
    } else {
        let temp =
            (julian_day as f64 - day_points[idx - 1]) / (day_points[idx] - day_points[idx - 1]);
        return irradiance_points[idx - 1]
            - (irradiance_points[idx - 1] - irradiance_points[idx]) * temp;
    }
}

fn compute_airmass(zen_r: f64, zen_d: f64) -> f64 {
    let mut airmass: f64;

    airmass = 1.0 / (zen_r.cos() + 0.15 * (93.885 - zen_d).powf(-1.253));

    if airmass < 1.0 {
        airmass = 1.0;
    }

    airmass
}

fn compute_rayleigh(airmass: f64) -> [f64; TRANSMITTANCE_WL_COUNT] {
    let mut tr: [f64; TRANSMITTANCE_WL_COUNT] = [0.0; TRANSMITTANCE_WL_COUNT];

    for w in 0..TRANSMITTANCE_WL_COUNT {
        let wld = TRANSMITTANCE_WAVELENGTHS[w] / 1000.0;
        tr[w] = (-airmass / (wld.powf(4.0) * (115.6406 - 1.335 / wld.powf(2.0)))).exp();
    }

    tr
}

fn compute_aerosol_transmittance(airmass: f64) -> [f64; TRANSMITTANCE_WL_COUNT] {
    let mut ta: [f64; TRANSMITTANCE_WL_COUNT] = [0.0; TRANSMITTANCE_WL_COUNT];

    for wl in 0..TRANSMITTANCE_WL_COUNT {
        let wld = TRANSMITTANCE_WAVELENGTHS[wl] / 1000.0; // NOTE: This is repeated a few times, worth a single pre-compute?
        if wl < 10 {
            ta[wl] = (-BETA1 * wld.powf(-ALPHA1) * airmass).exp();
        } else {
            ta[wl] = (-BETA2 * wld.powf(-ALPHA2) * airmass).exp();
        }
    }
    ta
}

fn compute_water_vapour_transmittance(airmass: f64) -> [f64; TRANSMITTANCE_WL_COUNT] {
    let mut tw: [f64; TRANSMITTANCE_WL_COUNT] = [0.0; TRANSMITTANCE_WL_COUNT];

    for w in 0..TRANSMITTANCE_WL_COUNT {
        tw[w] = (-0.3285 * WATER_VAPOUR_ABS[w] * (W + (1.42 - W) / 2.0) * airmass
            / (1.0 + 20.07 * WATER_VAPOUR_ABS[w] * airmass).powf(0.45))
        .exp()
    }
    tw
}

fn compute_ozone_transmittance(zenith_r: f64) -> [f64; TRANSMITTANCE_WL_COUNT] {
    let mut to: [f64; TRANSMITTANCE_WL_COUNT] = [0.0; TRANSMITTANCE_WL_COUNT];

    for w in 0..TRANSMITTANCE_WL_COUNT {
        let em0: f64 = 35.0 / (1224.0 * (zenith_r.cos()).powf(2.0) + 1.0).powf(0.5);
        to[w] = (-OZONE_ABS[w] * 0.03 * em0).exp();
    }

    to
}

fn compute_tu(airmass: f64) -> f64 {
    (-1.41 * 0.15 * airmass / (1.0 + 118.3 * 0.15 * airmass).powf(0.45)).exp()
}

fn compute_air_albedo(
    aerosol: [f64; TRANSMITTANCE_WL_COUNT],
    ozone: [f64; TRANSMITTANCE_WL_COUNT],
    rayleigh: [f64; TRANSMITTANCE_WL_COUNT],
    water_vapour: [f64; TRANSMITTANCE_WL_COUNT],
    tu: f64,
) -> [f64; TRANSMITTANCE_WL_COUNT] {
    let mut air_albedo: [f64; TRANSMITTANCE_WL_COUNT] = [0.0; TRANSMITTANCE_WL_COUNT];

    for w in 0..TRANSMITTANCE_WL_COUNT {
        air_albedo[w] = ozone[w]
            * water_vapour[w]
            * (aerosol[w] * (1.0 - rayleigh[w]) * 0.5
                + rayleigh[w] * (1.0 - aerosol[w]) * 0.22 * 0.928);
    }

    air_albedo[23] = air_albedo[23] * tu;
    air_albedo
}

// given the calculated zenith angle (DEGREES), look up the index of the correction
// factors we need to use
fn find_lut_index(zen_d: f64) -> usize {
    for i in 0..CORRECTION_ZEN_LOOKUP.len() {
        if zen_d < CORRECTION_ZEN_LOOKUP[i] {
            return i;
        }
    }
    return CORRECTION_ZEN_LOOKUP.len() - 1;
}

// take the correction factors for 5 wavelengths from the LUT and interpolate to the
// number of wavelengths we are using in the transmittance calculations (24 at time of writing)
fn interpolate_correction_factor(zen_d: f64) -> [f64; TRANSMITTANCE_WL_COUNT] {
    let lut_index: usize = find_lut_index(zen_d);

    let c: [f64; 5] = CORRECTION[lut_index];
    let cm1: [f64; 5] = CORRECTION[lut_index - 1];

    let fraction: f64 = (zen_d - CORRECTION_ZEN_LOOKUP[lut_index - 1])
        / (CORRECTION_ZEN_LOOKUP[lut_index] - CORRECTION_ZEN_LOOKUP[lut_index - 1]);

    let mut correction: [f64; 5] = [0.0; 5];
    for i in 0..c.len() {
        correction[i] = (c[i] - cm1[i]) * fraction + cm1[i];
    }
    let correction = correction;

    let mut correction_interpolated: [f64; TRANSMITTANCE_WL_COUNT] = [0.0; TRANSMITTANCE_WL_COUNT];

    let mut c_inc: f64;
    let mut l: usize = 0;

    correction_interpolated[0] = correction[0];

    for l1 in 1..4 {
        c_inc = (correction[l1 - 1] - correction[l1]) / 5.0;
        for _l2 in 0..5 {
            l += 1;
            correction_interpolated[l] = correction_interpolated[l - 1] - c_inc;
        }
    }

    // for the remaining wavelengths 550 -> 710nm
    let c_dif: f64 = correction[4] - correction[3];
    let wldif =
        TRANSMITTANCE_WAVELENGTHS[TRANSMITTANCE_WL_COUNT - 1] - TRANSMITTANCE_WAVELENGTHS[15];

    for l1 in 16..24 {
        let wl_inc = (TRANSMITTANCE_WAVELENGTHS[l1] - TRANSMITTANCE_WAVELENGTHS[l1 - 1]) / wldif;
        c_inc = c_dif * wl_inc;
        l += 1;
        correction_interpolated[l] = correction_interpolated[l - 1] + c_inc;
    }

    return correction_interpolated;
}

fn compute_diffuse_irradiance(
    zen_r: f64,
    zen_d: f64,
    direct_irradiance: [f64; TRANSMITTANCE_WL_COUNT],
    air_albedo: [f64; TRANSMITTANCE_WL_COUNT],
    t_aerosol: [f64; TRANSMITTANCE_WL_COUNT],
    t_ozone: [f64; TRANSMITTANCE_WL_COUNT],
    t_rayleigh: [f64; TRANSMITTANCE_WL_COUNT],
    tu: f64,
    t_water_vapour: [f64; TRANSMITTANCE_WL_COUNT],
) -> [f64; TRANSMITTANCE_WL_COUNT] {
    let correction_matrix: [f64; TRANSMITTANCE_WL_COUNT] = interpolate_correction_factor(zen_d);

    let mut diffuse: [f64; TRANSMITTANCE_WL_COUNT] = [0.0; TRANSMITTANCE_WL_COUNT];

    for l in 0..24 {
        let xx = ET_SPECTRAL_IRRADIANCE[l] * zen_r.cos() * t_ozone[l] * t_water_vapour[l];
        let mut r = xx * t_aerosol[l] * (1.0 - t_rayleigh[l]) * 0.5;
        let mut a = xx * t_rayleigh[l] * (1.0 - t_aerosol[l] * 0.928 * 0.82);

        if l == 22 {
            r = r * tu;
            a = a * tu;
        }

        let g = (direct_irradiance[l] * zen_r.cos() + (r + a) * correction_matrix[l])
            * air_albedo[l]
            * 0.05
            / (1.0 - 0.05 * air_albedo[l]);

        diffuse[l] = (r + a) * correction_matrix[l] + g;
    }

    return diffuse;
}

fn compute_direct_irradiance(
    t_aerosol: [f64; TRANSMITTANCE_WL_COUNT],
    t_ozone: [f64; TRANSMITTANCE_WL_COUNT],
    t_rayleigh: [f64; TRANSMITTANCE_WL_COUNT],
    tu: f64,
    t_water_vapour: [f64; TRANSMITTANCE_WL_COUNT],
) -> [f64; TRANSMITTANCE_WL_COUNT] {
    let mut direct: [f64; TRANSMITTANCE_WL_COUNT] = [0.0; TRANSMITTANCE_WL_COUNT];

    for i in 0..TRANSMITTANCE_WL_COUNT {
        direct[i] = ET_SPECTRAL_IRRADIANCE[i]
            * t_rayleigh[i]
            * t_aerosol[i]
            * t_water_vapour[i]
            * t_ozone[i]
    }
    direct[23] *= tu;

    direct
}

pub fn compute_irradiance_components(
    zenith_r: f64,
    zenith_d: f64,
) -> ([f64; WL_COUNT], [f64; WL_COUNT]) {
    // use airmass estimate initially until we calculate air albedo and then we recalculate transmittances
    let airmass = 1.90;

    let t_rayleigh = compute_rayleigh(airmass);
    let t_aerosol = compute_aerosol_transmittance(airmass);
    let t_water_vapour = compute_water_vapour_transmittance(airmass);
    let t_ozone = compute_ozone_transmittance(zenith_r);
    let t_u = compute_tu(airmass);

    let air_albedo = compute_air_albedo(t_aerosol, t_ozone, t_rayleigh, t_water_vapour, t_u);

    let airmass = compute_airmass(zenith_r, zenith_d);
    let t_rayleigh = compute_rayleigh(airmass);
    let t_aerosol = compute_aerosol_transmittance(airmass);
    let t_water_vapour = compute_water_vapour_transmittance(airmass);
    let t_ozone = compute_ozone_transmittance(zenith_r);
    let t_u = compute_tu(airmass);

    let direct = compute_direct_irradiance(t_aerosol, t_ozone, t_rayleigh, t_u, t_water_vapour);

    let diffuse = compute_diffuse_irradiance(
        zenith_r,
        zenith_d,
        direct,
        air_albedo,
        t_aerosol,
        t_ozone,
        t_rayleigh,
        t_u,
        t_water_vapour,
    );

    let direct_interpolated = interpolate_irradiances(TRANSMITTANCE_WAVELENGTHS, WL_ARRAY, direct);
    let diffuse_interpolated =
        interpolate_irradiances(TRANSMITTANCE_WAVELENGTHS, WL_ARRAY, diffuse);

    return (direct_interpolated, diffuse_interpolated);
}

fn interpolate_irradiances(
    input_wavelengths: [f64; TRANSMITTANCE_WL_COUNT],
    output_wavelengths: [f64; WL_COUNT],
    input_irradiances: [f64; TRANSMITTANCE_WL_COUNT],
) -> [f64; WL_COUNT] {
    let mut output_irradiances: [f64; WL_COUNT] = [0.0; WL_COUNT];

    for i in 0..WL_COUNT {
        output_irradiances[i] = linear_interp(
            &input_wavelengths,
            &input_irradiances,
            output_wavelengths[i],
        );
    }

    output_irradiances
}
