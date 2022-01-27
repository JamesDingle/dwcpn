pub mod dwcpn;

#[cfg(test)]
mod tests {
    use crate::dwcpn::dwcpn::{calc_pp, ModelInputs, ModelSettings};
    use crate::dwcpn::modules::pp_profile::{calculate_ay, calculate_bbr, calculate_bw};
    
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
        
        let (pp, _, _) = calc_pp(&inputs, &settings);
        println!("PP Test result: {} (expected 721.7......)", pp);
        assert!(pp >= 721.0);
        assert!(pp <= 722.0);
    }
}
