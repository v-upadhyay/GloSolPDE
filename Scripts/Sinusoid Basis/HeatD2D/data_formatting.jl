using MAT
using DelimitedFiles
using Serialization

function convert_dat_to_mat(dat_file::String, mat_file::String)
    # Read the .dat file
    data = deserialize(open(dat_file, "r"))
    
    # Save the data to a .mat file
    matwrite(mat_file, Dict("data" => data))
end

# Example usage
convert_dat_to_mat("recent_traj.dat", "Trajectory.mat")
