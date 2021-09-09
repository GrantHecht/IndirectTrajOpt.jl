
# Function for reading binary data files (.jld2) into memory. Returns vector
# of IndirectTrajOptimizer(s).
#
# folder:   String containing the full file path to folder containing binary
#           data files.
#

function readBinaryData(folder::String)
    # Get all files in folder
    files = readdir(joinpath(folder, "binaryData"))

    # Read all files in push data to data vector
    local dataVec
    @showprogress "Reading Binary Data Files: " for i in 1:length(files)
        # Open file
        f = jldopen(joinpath(folder, "binaryData", files[i]), "r")

        # Grab data and push to vector
        if i == 1
            dataVec = [f["data"]]
        else
            push!(dataVec, f["data"])
        end

        # Close file
        close(f)
    end

    return dataVec
end