use crate::dwcpn::modules::config::{DEPTH_PROFILE_COUNT, WL_COUNT};

pub mod dwcpn;

pub struct ModelInputs {
    pub lat: f64,
    pub lon: f64,
    pub z_bottom: f64,
    pub iday: u16,
    pub alpha_b: f64,
    pub pmb: f64,
    pub z_m: f64,
    pub mld: f64,
    pub chl: f64,
    pub rho: f64,
    pub sigma: f64,
    pub cloud: f64,
    pub yel_sub: f64,
    pub par: f64,
    pub bw: [f64; WL_COUNT],
    pub bbr: [f64; WL_COUNT],
    pub ay: [f64; WL_COUNT],
}

pub struct ProchloroInputs {
    pub prochloro_surface: f64,
    pub prochloro_maximum: f64
}

pub struct ModelSettings {
    pub mld_only: bool,
    pub iom_only: bool,
    pub prochloro_inputs: Option<ProchloroInputs>
}

pub struct ModelOutputs {
    pub pp_day: Option<f64>,
    pub euphotic_depth: Option<f64>,
    pub spectral_i_star: Option<f64>,
    pub par_noon_max: Option<f64>,
    pub pro_1_profile: Option<[f64; DEPTH_PROFILE_COUNT]>,
    pub pro_2_profile: Option<[f64; DEPTH_PROFILE_COUNT]>,
    pub pro_total_profile: Option<[f64; DEPTH_PROFILE_COUNT]>,
}

#[derive(Debug)]
pub enum PPErrors {
    DWCPNError,
    PPTooBig
}

#[cfg(test)]
mod tests {
    use crate::dwcpn::dwcpn::{calc_production};
    use crate::dwcpn::modules::pp_profile::{calculate_ay, calculate_bbr, calculate_bw, compute_prochloro_profile};
    use crate::{ModelInputs, ModelSettings, ProchloroInputs};

    #[test]
    fn pp_integration_test() {
        
        let settings = ModelSettings {
            mld_only: false,
            iom_only: false,
            prochloro_inputs: None
        };
        
        let inputs = ModelInputs {
            lat: -5.792,
            lon: -96.62,
            z_bottom: 100.0,
            iday: 1,
            alpha_b: 0.0844,
            pmb: 4.756,
            z_m: 46.1,
            mld: 19.35091019,
            chl: 0.26096588,
            rho: 0.878,
            sigma: 34.6,
            cloud: 0.0,
            yel_sub: 0.3,
            par: 49.1697464,
            bw: calculate_bw(),
            bbr: calculate_bbr(),
            ay: calculate_ay()
        };

        match calc_production(&inputs, &settings) {
            Ok(model_output) => {
                println!("PP Test result: {} (expected 721.7......)", model_output.pp_day.unwrap());
                println!("IOM output: {}", model_output.par_noon_max.unwrap());
                assert!(model_output.pp_day.unwrap() >= 721.0);
                assert!(model_output.pp_day.unwrap() <= 722.0);
            },
            Err(e) => {
                println!("{:?}", e);
            }
        }

    }

    #[test]
    fn prochloro_integration_test() {

        let settings = ModelSettings {
            mld_only: false,
            iom_only: false,
            prochloro_inputs: Some(ProchloroInputs {
                prochloro_surface: 1.0,
                prochloro_maximum: 2.0
            })
        };

        let inputs = ModelInputs {
            lat: -5.792,
            lon: -96.62,
            z_bottom: 100.0,
            iday: 1,
            alpha_b: 0.0844,
            pmb: 4.756,
            z_m: 46.1,
            mld: 19.35091019,
            chl: 0.26096588,
            rho: 0.878,
            sigma: 34.6,
            cloud: 0.0,
            yel_sub: 0.3,
            par: 49.1697464,
            bw: calculate_bw(),
            bbr: calculate_bbr(),
            ay: calculate_ay()
        };

        match calc_production(&inputs, &settings) {
            Ok(model_output) => {
                println!("PP Test result: {} (expected 721.7......)", model_output.pp_day.unwrap());
                println!("IOM output: {}", model_output.par_noon_max.unwrap());
                assert!(model_output.pp_day.unwrap() >= 721.0);
                assert!(model_output.pp_day.unwrap() <= 722.0);
            },
            Err(e) => {
                println!("{:?}", e);
            }
        }

    }
}
