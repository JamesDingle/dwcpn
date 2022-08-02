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
    pub pp_prochloro_profile: Option<[f64; DEPTH_PROFILE_COUNT]>
}

#[derive(Debug)]
pub enum PPErrors {
    DWCPNError,
    PPTooHigh
}

#[cfg(test)]
mod integration_tests {
    use crate::dwcpn::dwcpn::{calc_production};
    use crate::dwcpn::modules::pp_profile::{calculate_ay, calculate_bbr, calculate_bw};
    use crate::{ModelInputs, ModelOutputs, ModelSettings};

    struct TestCase {
        name: String,
        inputs: ModelInputs,
        settings: ModelSettings,
        expected_result: ModelOutputs
    }

    fn run_test_case(test_case: &TestCase, accuracy: f64) -> bool {
        match calc_production(&test_case.inputs, &test_case.settings) {
            Ok(model_output) => {

                let result = model_output.pp_day.unwrap();
                let expected_result = test_case.expected_result.pp_day.unwrap();

                let lower_bound = expected_result / 100.0 * (100.0 - accuracy);
                let upper_bound = expected_result / 100.0 * (100.0 + accuracy);

                if result >= lower_bound && result <= upper_bound {
                    println!("{}: Success! Result {:.3}, (target:{:.2}<-{:.2}->{:.2})",
                             test_case.name,
                             result,
                             lower_bound,
                             expected_result,
                             upper_bound
                    );
                    true
                } else {
                    println!("{}: Fail! Result {:.3}, (target:{:.2}<-{:.2}->{:.2})",
                             test_case.name,
                             result,
                             lower_bound,
                             expected_result,
                             upper_bound
                    );
                    false
                }
            },
            Err(e) => {
                println!("{:?}", e);
                false
            }
        }
    }

    #[test]
    fn east_pacific_test_pp() {
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

        let settings = ModelSettings {
            mld_only: false,
            iom_only: false,
            prochloro_inputs: None
        };

        let expected_result = ModelOutputs {
            pp_day: Some(721.7),
            euphotic_depth: None,
            spectral_i_star: None,
            par_noon_max: None,
            pro_1_profile: None,
            pro_2_profile: None,
            pro_total_profile: None,
            pp_prochloro_profile: None
        };

        let test_case = TestCase {
            name: "east_pacific".to_string(),
            inputs,
            settings,
            expected_result
        };

        assert!(run_test_case(&test_case, 2.0));

    }

    #[test]
    fn black_sea_test_pp() {
        let inputs = ModelInputs {
            lat: 43.2,
            lon: 33.7,
            z_bottom: 2198.8,
            iday: 121,
            alpha_b: 0.0578,
            pmb: 3.294,
            z_m: 49.44,
            mld: 11.9296,
            chl: 0.474,
            rho: 0.87,
            sigma: 14.62,
            cloud: 0.0,
            yel_sub: 0.3,
            par: 50.35,
            bw: calculate_bw(),
            bbr: calculate_bbr(),
            ay: calculate_ay()
        };

        let settings = ModelSettings {
            mld_only: false,
            iom_only: false,
            prochloro_inputs: None
        };

        let expected_result = ModelOutputs {
            pp_day: Some(905.18976),
            euphotic_depth: None,
            spectral_i_star: None,
            par_noon_max: None,
            pro_1_profile: None,
            pro_2_profile: None,
            pro_total_profile: None,
            pp_prochloro_profile: None
        };

        let test_case = TestCase {
            name: "black_sea".to_string(),
            inputs,
            settings,
            expected_result
        };

        assert!(run_test_case(&test_case, 5.0));
    }

    #[test]
    fn south_atlantic_gyre_test_pp() {
        let inputs = ModelInputs {
            lat: -27.042,
            lon: -17.71,
            z_bottom: 3940.65,
            iday: 121,
            alpha_b: 0.0933,
            pmb: 1.594,
            z_m: 97.82,
            mld: 62.349,
            chl: 0.058,
            rho: 0.75,
            sigma: 20.64,
            cloud: 0.0,
            yel_sub: 0.3,
            par: 25.482,
            bw: calculate_bw(),
            bbr: calculate_bbr(),
            ay: calculate_ay()
        };

        let settings = ModelSettings {
            mld_only: false,
            iom_only: false,
            prochloro_inputs: None
        };

        let expected_result = ModelOutputs {
            pp_day: Some(108.63),
            euphotic_depth: None,
            spectral_i_star: None,
            par_noon_max: None,
            pro_1_profile: None,
            pro_2_profile: None,
            pro_total_profile: None,
            pp_prochloro_profile: None
        };

        let test_case = TestCase {
            name: "south_atlantic_gyre".to_string(),
            inputs,
            settings,
            expected_result
        };

        assert!(run_test_case(&test_case, 5.0));
    }

    #[test]
    fn mauritania_upwelling_test_pp() {
        let inputs = ModelInputs {
            lat: 18.71,
            lon: -18.625,
            z_bottom: 2950.468,
            iday: 121,
            alpha_b: 0.1518,
            pmb: 3.9059,
            z_m: 23.094,
            mld: 31.975,
            chl: 1.718,
            rho: 0.8247,
            sigma: 27.556,
            cloud: 0.0,
            yel_sub: 0.3,
            par: 55.8677,
            bw: calculate_bw(),
            bbr: calculate_bbr(),
            ay: calculate_ay()
        };

        let settings = ModelSettings {
            mld_only: false,
            iom_only: false,
            prochloro_inputs: None
        };

        let expected_result = ModelOutputs {
            pp_day: Some(2341.988),
            euphotic_depth: None,
            spectral_i_star: None,
            par_noon_max: None,
            pro_1_profile: None,
            pro_2_profile: None,
            pro_total_profile: None,
            pp_prochloro_profile: None
        };

        let test_case = TestCase {
            name: "mauritania upwelling".to_string(),
            inputs,
            settings,
            expected_result
        };

        assert!(run_test_case(&test_case, 1.0));
    }

    #[test]
    fn arabian_sea_test_pp() {
        let inputs = ModelInputs {
            lat: 12.542,
            lon: 62.542,
            z_bottom: 4244.75,
            iday: 121,
            alpha_b: 0.1329,
            pmb: 3.952,
            z_m: 70.42,
            mld: 26.15,
            chl: 0.1032,
            rho: 0.856,
            sigma: 23.523,
            cloud: 0.0,
            yel_sub: 0.3,
            par: 56.255,
            bw: calculate_bw(),
            bbr: calculate_bbr(),
            ay: calculate_ay()
        };

        let settings = ModelSettings {
            mld_only: false,
            iom_only: false,
            prochloro_inputs: None
        };

        let expected_result = ModelOutputs {
            pp_day: Some(694.43),
            euphotic_depth: None,
            spectral_i_star: None,
            par_noon_max: None,
            pro_1_profile: None,
            pro_2_profile: None,
            pro_total_profile: None,
            pp_prochloro_profile: None
        };

        let test_case = TestCase {
            name: "arabian sea".to_string(),
            inputs,
            settings,
            expected_result
        };

        assert!(run_test_case(&test_case, 0.1));
    }

}
