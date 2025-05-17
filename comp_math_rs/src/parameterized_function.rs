pub enum ParameterizedFunctionError {
    WrongParameterId,
}

pub trait ParameterizedFunction<const ARGS_N: usize, const RES_N: usize> {
    fn calc(&mut self, args: &[f64; ARGS_N]) -> [f64; RES_N];
    fn set_parameter(&mut self, id: usize, parameter: f64) -> Result<(), ParameterizedFunctionError>;
    fn get_parameter(&self, id: usize) -> Result<f64, ParameterizedFunctionError>;
}
