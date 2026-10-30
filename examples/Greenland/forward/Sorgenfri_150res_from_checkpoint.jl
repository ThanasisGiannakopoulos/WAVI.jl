using WAVI
using FileIO

function main()

    chkpt_file = "./PChkpt_0000001600.jld2"

    chk = load(chkpt_file)

    println("Checkpoint contents:")
    println(keys(chk))

    # Replace with the actual key present in your checkpoint
    simulation = chk["simulation"]

    run_simulation!(simulation)

end

main()