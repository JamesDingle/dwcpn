pub mod dwcpn;

#[cfg(test)]
mod tests {
    use crate::dwcpn::dwcpn::{calc_pp, ModelInputs, ModelSettings};
    use crate::dwcpn::modules::pp_profile::{calculate_ay, calculate_bbr, calculate_bw};

    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
    
    #[test]
    fn test() {
        
        let settings = ModelSettings {
            mld_only: false
        };
        
        let inputs = ModelInputs {
            lat: 50.,
            lon: 0.0,
            z_bottom: 100.0,
            iday: 1,
            alpha_b: 1.0,
            pmb: 1.0,
            z_m: 1.0,
            mld: 20.0,
            chl: 1.0,
            rho: 1.0,
            sigma: 1.0,
            cloud: 0.0,
            yel_sub: 0.3,
            par: 100.0,
            bw: calculate_bw(),
            bbr: calculate_bbr(),
            ay: calculate_ay()
        };
        
        let (pp, _, _) = calc_pp(&inputs, &settings);
        println!("Test pp result: {}", pp);
    }
}
