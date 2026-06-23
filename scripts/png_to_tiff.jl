#!/usr/bin/env julia

using FileIO

function print_usage()
    println("Usage:")
    println("  julia scripts/png_to_tiff.jl <folder> [--recursive] [--overwrite]")
    println()
    println("Options:")
    println("  --recursive   Search for .png files in subfolders as well")
    println("  --overwrite   Overwrite existing .tiff files")
end

function is_png(path::AbstractString)
    return lowercase(last(splitext(path))) == ".png"
end

function find_pngs(folder::AbstractString; recursive::Bool=false)
    pngs = String[]

    if recursive
        for (root, _, files) in walkdir(folder)
            for file in files
                fullpath = joinpath(root, file)
                is_png(fullpath) && push!(pngs, fullpath)
            end
        end
    else
        for name in readdir(folder)
            fullpath = joinpath(folder, name)
            if isfile(fullpath) && is_png(fullpath)
                push!(pngs, fullpath)
            end
        end
    end

    return sort(pngs)
end

function convert_pngs(folder::AbstractString; recursive::Bool=false, overwrite::Bool=false)
    if !isdir(folder)
        error("Folder does not exist: $folder")
    end

    png_files = find_pngs(folder; recursive)

    if isempty(png_files)
        println("No .png files found in: $folder")
        return
    end

    converted = 0
    skipped = 0
    failed = 0

    for png in png_files
        base, _ = splitext(png)
        tiff = base * ".tiff"

        if isfile(tiff) && !overwrite
            println("[skip] $tiff already exists")
            skipped += 1
            continue
        end

        try
            img = load(png)
            save(tiff, img)
            println("[ok]   $png -> $tiff")
            converted += 1
        catch err
            println("[fail] $png")
            println("       $(sprint(showerror, err))")
            failed += 1
        end
    end

    println()
    println("Done.")
    println("Converted: $converted")
    println("Skipped:   $skipped")
    println("Failed:    $failed")
end

function main(args)
    if isempty(args) || "--help" in args || "-h" in args
        print_usage()
        return
    end

    folder = args[1]
    recursive = "--recursive" in args
    overwrite = "--overwrite" in args

    convert_pngs(folder; recursive, overwrite)
end

main(ARGS)
