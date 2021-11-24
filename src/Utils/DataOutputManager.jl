
mutable struct DataOutputManager

    # Base folder location. Should be full file path and not local
    baseFolder::String

    # File name without extension
    fname::String

    # File name initialized flag
    fnameInitialized::Bool

    # Text file format version
    textFormatVersion::String

    # Constructor
    function DataOutputManager(baseFolder::String)
        new(baseFolder, "", false, "v0.1")
    end
end

