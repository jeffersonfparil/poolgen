dir = "/home/jeffersonfparil/Documents/poolgen/test"
cd(dir)

vec_files = readdir()
for f in vec_files[match.(Regex("FILTERED"), vec_files) .!= nothing]
    rm(f)
end

vec_files = readdir()
for f in vec_files[match.(Regex("IMPUTED"), vec_files) .!= nothing]
    rm(f)
end

vec_files = readdir()
for f in vec_files[match.(Regex("RAW"), vec_files) .!= nothing]
    rm(f)
end

vec_files = readdir()
for f in vec_files[match.(Regex("syncx"), vec_files) .!= nothing]
    rm(f)
end
