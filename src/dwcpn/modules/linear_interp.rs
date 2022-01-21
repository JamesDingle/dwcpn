// I have put the interpolation into its own module so that I can experiment with different methods
// in future

pub fn linear_interp(index_array: &[f64], value_array: &[f64], index_lookup: f64) -> f64 {
    let size = index_array.len();

    let mut i: usize = 0;

    if index_lookup >= index_array[size - 2] {
        i = size - 2;
    } else {
        while index_lookup > index_array[i + 1] {
            i = i + 1
        }
    }

    let index_left = index_array[i];
    let index_right = index_array[i + 1];

    let value_left = value_array[i];
    let value_right = value_array[i + 1];

    if (index_lookup < index_left) || (index_lookup > index_right) {
        println!(
            "{:?} < {:?} or {:?} > {:?}",
            index_lookup, index_left, index_lookup, index_right
        );
        panic!("Trying to interpolate outside of index array. Extrapolation is forbidden in this use case");
    }

    let gradient = (value_right - value_left) / (index_right - index_left);

    let result = value_left + gradient * (index_lookup - index_left);

    return result;
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_non_exact() {
        let x: [f64; 5] = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y: [f64; 5] = [2.0, 4.0, 6.0, 8.0, 10.0];
        let lookup: f64 = 2.5;

        assert_eq!(linear_interp(&x, &y, lookup), 5.0);
    }

    #[test]
    fn test_exact() {
        let x: [f64; 5] = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y: [f64; 5] = [2.0, 4.0, 6.0, 8.0, 10.0];
        let lookup: f64 = 3.0;

        assert_eq!(linear_interp(&x, &y, lookup), 6.0);
    }

    #[test]
    fn test_exact_first_element() {
        let x: [f64; 5] = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y: [f64; 5] = [2.0, 4.0, 6.0, 8.0, 10.0];

        let lookup: f64 = 1.0;
        assert_eq!(linear_interp(&x, &y, lookup), 2.0);
    }

    #[test]
    fn test_exact_last_element() {
        let x: [f64; 5] = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y: [f64; 5] = [2.0, 4.0, 6.0, 8.0, 10.0];

        let lookup: f64 = 5.0;
        assert_eq!(linear_interp(&x, &y, lookup), 10.0);
    }

    #[test]
    #[should_panic]
    fn test_extrapolation() {
        let x: [f64; 5] = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y: [f64; 5] = [2.0, 4.0, 6.0, 8.0, 10.0];

        let lookup: f64 = 6.0;
        linear_interp(&x, &y, lookup);
    }
}
