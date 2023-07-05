"""
    Make folder that have non-repeating index by the end of its name. Give the folder an "sID" index (simulation Index) automatically if not specified.

    Returns the full path to the folder being created.
"""
function make_indexed_folder(; folder_prefix, folder_path=".", sID::Int=0)

    # initialize data folder names
    folder_name = @sprintf "%s-%d" folder_prefix sID
    folder_full_path = joinpath(folder_path, folder_name)

    # if null data folder id given, determine data name and id
    if sID == 0
        while isdir(folder_full_path) || sID == 0
            sID += 1
            folder_name = @sprintf "%s-%d" folder_prefix sID
            folder_full_path = joinpath(folder_path, folder_name)
        end
    end

    # make the folder
    mkpath(folder_full_path)

    return folder_full_path
end
