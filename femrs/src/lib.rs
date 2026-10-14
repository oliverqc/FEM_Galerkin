struct InformationMatrices {
    pub elements: u16, 
    pub dimension: u16, 
}

struct BasisFunctions {
    pub basis_type: String, 
    pub num_functions: u16,
}

impl InformationMatrices {
    // Builds T information matrix for linear, quadratic, 1D, 2D.
    fn mesh_elements(&self, basis_functions: &BasisFunctions) -> Option<Vec<u16>> {
        let vec: Vec<f32> = Vec::new(); 
        if basis_functions.num_functions == 2 {
            vec.extend(1..self.elements-1);
            vec.extend(2..self.elements);
            Some(vec)
        } else if {}
    }
    // Creates P information matrix. 
    fn mesh_coordinates(&self) -> Option<Vec<f32>> {
        todo!();
    }
}
