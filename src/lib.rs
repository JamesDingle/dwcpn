use crate::dwcpn::modules::config::WL_COUNT;

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

pub struct ModelSettings {
    pub mld_only: bool
}

pub struct ModelOutputs {
    pub pp_day: Option<f64>,
    pub euphotic_depth: Option<f64>,
    pub spectral_i_star: Option<f64>,
    pub par_noon_max: Option<f64>
}

#[derive(Debug)]
pub enum PPErrors {
    DWCPNError,
    PPTooBig
}

#[cfg(test)]
mod tests {
    use crate::dwcpn::dwcpn::{calc_pp};
    use crate::dwcpn::modules::pp_profile::{calculate_ay, calculate_bbr, calculate_bw};
    use crate::{ModelInputs, ModelSettings};

    #[test]
    fn integration_test() {
        
        let settings = ModelSettings {
            mld_only: false
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

        match calc_pp(&inputs, &settings) {
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
