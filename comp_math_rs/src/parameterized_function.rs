pub enum ParameterizedFunctionError {
    WrongParameterId,
}

pub trait ParameterizedFunction<const ARGS_N: usize, const RES_N: usize> {
    fn calc(&mut self, args: &[f64; ARGS_N]) -> [f64; RES_N];
    fn set_parameter(
        &mut self,
        id: usize,
        parameter: f64,
    ) -> Result<(), ParameterizedFunctionError>;
    fn get_parameter(&self, id: usize) -> Result<f64, ParameterizedFunctionError>;
}

#[macro_export]
macro_rules! parameterized_function {
    ($input_dim:literal, $output_dim:literal, $name:ident, $body:expr) => {
        struct $name {}

        impl $crate::ParameterizedFunction<$input_dim, $output_dim> for $name {
            fn calc(&mut self, x: &[f64; $input_dim]) -> [f64; $output_dim] {
                $body
            }

            fn get_parameter(&self, _id: usize) -> Result<f64, $crate::ParameterizedFunctionError> {
                Err($crate::ParameterizedFunctionError::WrongParameterId)
            }

            fn set_parameter(
                &mut self,
                _id: usize,
                _parameter: f64,
            ) -> Result<(), $crate::ParameterizedFunctionError> {
                Err($crate::ParameterizedFunctionError::WrongParameterId)
            }
        }
    };
    ($input_dim:literal, $output_dim:literal, $body:expr) => {
        parameterized_function!($input_dim, $output_dim, AnonymousFunction, $body);
    };
}
