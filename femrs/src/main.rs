use femrs::{InformationMatrices, BasisFunctions};

fn main() {
    let information_matrice = InformationMatrices {
        elements: 11,
        dimension: 1,
    };

    let linear_basis = BasisFunctions {
        num_functions: 2,
    };

    let mesh_element_vec = information_matrice.mesh_elements(&linear_basis).unwrap();
    println!("{:?}", mesh_element_vec);

    let mesh_coordinate_vec = information_matrice.mesh_coordinates().unwrap(); 
    println!("{:?}", mesh_coordinate_vec);
}





