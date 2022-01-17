use std::io::{BufRead, BufReader};
use std::fs;
use std::fs::File;
use std::env;
use std::io::{Write};
use std::str::FromStr;
fn load_data(filename: &str) -> std::io::Result<Vec<Vec<f64>>> {
    let mut data = Vec::new();
    let mut l = 0;

    for result in BufReader::new(fs::File::open(filename)?).lines() {
        let line = Vec::new();
        data.push(line);
        for s in result?.split(',') {
            // println!("{}", s);
            let l_s = match f64:: from_str(s) {
                Ok(value) => value,
                Err(_err) => 0.0,
            };
            data[l].push(l_s);
        }
        l += 1;
    }
    Ok(data)
}

fn main() {
    const pi         : f64 = std::f64::consts::PI; // 3.141592...[無次元]
    let   m_N_c2     : f64 = 939.1 * 1.602e-13; // 939.1 [MeV] *  1.602e-13 [J / MeV] = [J] = [kg m^2 s^-2]
    let   m_N        : f64 = 1.6715e-27;           // 1.7e-27 [kg]
    let   k_B        : f64 = 1.380649e-23;         // 1.4e-23 [J/K]
    let   hbar       : f64 = 1.054571e-34;         // [J s]
    let   hbar_c     : f64 = hbar * 2.998e8;         // [J m]
    let   mu         : f64 = 1.660539e-27;           // [kg]
    let   NA         : f64 = 6.02214076e26;          // [kg^-1]


    let A            : f64 = k_B.powf(3.) * 7.0 * pi * pi  * m_N / (hbar_c.powf(3.) * 45.);
    let B            : f64 = (3. * pi.powf(2.)).powf(2./3.) * k_B  * m_N.powf(1./3.)/ (3. * hbar_c);


    let Ye : Vec<String> = ["001", "003", "006", "009", "010", "013", "016", "019", "020", "023", "026", "029", "030"].iter().map(|s| s.to_string()).collect();
    let T  : Vec<String> = ["4e9", "7e9", "1e10", "4e10", "7e10", "1e11", "4e11", "7e11", "1e12"].iter().map(|s| s.to_string()).collect();
    let rho: Vec<String> = ["1e10", "4e10", "7e10", "1e11", "4e11", "7e11", "1e12", "4e12", "7e12", "1e13"].iter().map(|s| s.to_string()).collect();
    
    for Ye_str in &Ye {
        for T_str in &T {
            for rho_str in &rho {
                
                let mut Ye_cal  : f64 = (Ye_str.parse().unwrap()); Ye_cal /= 100.;
                let mut T_cal   : f64 = T_str.parse().unwrap();
                let mut rho_cal : f64 = rho_str.parse().unwrap(); rho_cal *= 1000.;
                
                let S_rad       : f64 = A * T_cal.powf(3.) / rho_cal;
                let S_dec       : f64 = B * Ye_cal.powf(2./3.) * T_cal / rho_cal.powf(1./3.);
                
                let input_data = "../../csv/NSE_Ye_".to_string() + Ye_str + "_T0_" + T_str + "_rho0_" + rho_str + ".csv";

                if let Ok(data) = load_data(&input_data) {
                    let data = load_data(&input_data).unwrap();
                    let mut S_ideal: f64 = 0.;
                    for nuclei in data {
                        // itr : 0            1 2 3 4    5       6  7  
                        // data: NaN(Element) Z A N mass solarMF MF IMF
                        let Xi = nuclei[6] as f64;
                        let mi = nuclei[4] as f64 * mu;
                        let mi_c2 = mi * 9e16;
                        let Yi = Xi / (mi * NA);
                        let ni = Xi * rho_cal / mi;
                        let LOG_NAKAMI = (1./ni) * ((mi_c2 * k_B * T_cal) / (2. * pi * hbar_c.powf(2.))).powf(1.5);
                        // println!("Xi : {:E}, mi : {:E}[kg], Yi : {:E}, ni : {:E}", Xi, mi, Yi, ni);
                        
                        // println!("{}", ((mu / (rho_value * 1000. * Xi)) * ((m_N_natural * 1.6 * 1e-13 * k_B * T_value) / (2. * pi * (hbar_c).powf(2.)).powf(1.5))).log10());
                        let Sci = 
                            if Xi == 0.0 || Xi < 1e-300 {0.0} 
                            else {Yi * (5./2. + LOG_NAKAMI.log10())};
                        S_ideal += Sci;
                    }
                    let S = S_rad + S_dec + S_ideal;
                    println!("Ye = {:.2}, T = {:E}, rho = {:E};    S_rad = {:.3E}, S_deg = {:.3E}, S_ideal = {:.3E},    S = {:.3E}", Ye_cal, T_cal, rho_cal, S_rad, S_dec, S_ideal, S);
                }

            }
        }
    }
}
