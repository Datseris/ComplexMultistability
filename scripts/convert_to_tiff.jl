# Because PLOS Climate is not smart enough to handle .png figures.
using DrWatson
@quickactivate "ComplexMultistability"
using FileIO, ImageMagick

folder = papersdir("figures")
for f in readdir(folder)
    inpath = joinpath(folder, f)
    (isfile(inpath) && endswith(lowercase(f), ".png")) || continue
    endswith(f, "_inv.png") && continue # skip inverted images
    outpath = joinpath(folder, splitext(f)[1] * ".tif")
    save(outpath, load(inpath))
end
