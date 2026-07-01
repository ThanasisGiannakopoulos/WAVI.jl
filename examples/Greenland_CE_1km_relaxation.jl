using WAVI, MAT, FileIO, Printf

function main()

    # Grid and boundary conditions
    nx = 580
    ny = 719
    nσ = 26 
    x0 = 238000.0
    y0 = -2.5586e6
    dx = 1050.0
    dy = 1050.0

    # Paths: static fields from original setup, evolved fields from inversion output
    static_basepath = "/home/tgia02/repos/WAVIConstructor.jl/examples/outputs/Greenland_CE_1km/"
    inversion_output_path = "/home/tgia02/repos/WAVI.jl/example_drivers/Rosie_examples/Greenland_CE_1km_βpower0.1_Bpower_shelf0.1_Bpower_grounded0.01/"
    inversion_iteration = 25

    # ============================================================================
    # Load static fields (these don't change during inversion)
    # ============================================================================
    hm=Array{Float64}(undef,nx,ny);
    read!(joinpath(static_basepath, "h_mask.bin"),hm);
    hm = map.(Bool, round.(Int, hm));

    u_iszero=Array{Float64}(undef,nx+1,ny);
    read!(joinpath(static_basepath, "u_iszero.bin"),u_iszero);
    u_iszero = map.(Bool, round.(Int, u_iszero));

    v_iszero=Array{Float64}(undef,nx,ny+1);
    read!(joinpath(static_basepath,"v_iszero.bin"),v_iszero);
    v_iszero = map.(Bool, round.(Int, v_iszero));

    sigma_grid=Array{Float64}(undef,nσ);
    read!(joinpath(static_basepath,"sigma_grid.bin"),sigma_grid);

    bed=Array{Float64}(undef,nx,ny);
    read!(joinpath(static_basepath, "bed.bin"),bed);

    temp=Array{Float64}(undef,nx,ny,nσ);
    read!(joinpath(static_basepath, "temps.bin"),temp);

    # ============================================================================
    # Create grid (after loading masks and sigma_grid)
    # ============================================================================
    grid = Grid(nx = nx,
                ny = ny,
                nσ = nσ,
                x0 = x0,
                y0 = y0,
                dx = dx,
                dy = dy,
                h_mask = hm,
                u_iszero = u_iszero,
                v_iszero = v_iszero,
                σ = sigma_grid)

    # ============================================================================
    # Load evolved fields from inversion .mat output
    # ============================================================================

    # Construct filename for the inversion iteration (format: outfile0000000025.mat)
    mat_filename = @sprintf("outfile%010d.mat", inversion_iteration)
    mat_filepath = joinpath(inversion_output_path, mat_filename)

    if !isfile(mat_filepath)
        error("Inversion output file not found: $mat_filepath")
    end

    println("Loading inversion output from: $mat_filepath")
    mat_data = matread(mat_filepath) #load(mat_filepath)  # Loads a dictionary of arrays

    # Extract evolved fields from .mat file
    # Note: Field names in .mat file match the output names from your inversion driver
    h = mat_data["h"]  # Thickness (nx, ny)
    u = mat_data["u"]  # U velocity (nx+1, ny) 
    v = mat_data["v"]  # V velocity (nx, ny+1)
    beta = mat_data["beta"]  # Beta field (nx, ny)
    grounded_fraction = mat_data["grounded_fraction"];  # Grounded fraction (nx, ny)

    # load damage and viscocity from checkpoint .jld2
    checkpoint_file = load(inversion_output_path*"PChkpt_0000000025.jld2")
    fields_3d = checkpoint_file["inversion_simulation"].model.fields.g3d

    damage = fields_3d.Φ
    viscosity = fields_3d.η

    # Default WAVI parameter values (if not specified in your inversion)
    density_ice = 918.0  # Default from WAVI
    density_ocean = 1028.0  # Default from WAVI
    sea_level_wrt_geoid = 0.0  # Default from WAVI

    bed_speed = checkpoint_file["inversion_simulation"].model.fields.gh.bed_speed#mat_data["bed_speed"]
    haf = checkpoint_file["inversion_simulation"].model.fields.gh.haf #mat_data["haf"]
    println("Loaded bed_speed and haf from .jld2 file")

    # Create mask (same as h_mask)
    # mask is a boolean array (nx, ny) where true = ice exists, false = no ice
    mask = hm

    # ============================================================================
    # Build model
    # ============================================================================

    #Initial conditions from inversion output

    initial_conditions = InitialConditions(
        initial_thickness = h,
        initial_viscosity = viscosity,
        initial_temperature = temp,
        initial_damage = damage,
        initial_u_veloc = u,
        initial_v_veloc = v
    )


    # Accumulation rate (static, from original setup)
    accumulation_rate=Array{Float64}(undef,nx,ny);
    read!(joinpath(static_basepath, "accumulation_data.bin"),accumulation_rate);
    #accumulation_rate.=ntoh.(accumulation_rate);

    # For relaxation runs, we do not evolve shelves
    evolveShelves = false

    params = Params(accumulation_rate = accumulation_rate,
                    #weertman_c = weertman_c,
                    evolveShelves = evolveShelves)

    # Use QuadraticMeltRate (update γT value as needed)
    #melt_rate = QuadraticMeltRate(γT = 5.3e-5)  # UPDATE THIS: set your melt rate parameter

    #solver parameters
    solver_params = SolverParams(maxiter_picard = 2)

    #parallel_spec = SharedMemorySpec(ngridsx = 2, ngridsy=4, overlap=1, niterations=1)

    model = Model(grid = grid, 
                bed_elevation = bed, 
                initial_conditions = initial_conditions, 
                params=params, 
                #parallel_spec = parallel_spec, 
                solver_params = solver_params,
                fracture = DruckerPragerPhaseField(),
                #   sliding_law = WeertmanSlidingLaw(),
                #   basal_hydrology = NoHydrology(),
                #   thermo_dynamics = NoThermoDynamics(),
                #   verbose = true
                )
                #melt_rate = melt_rate  # Commented out - no melt for relaxation

    #
    #timestepping parameters
    niter0 = 0
    step_thickness = false
    dt = 0.1
    end_time = 1. #20.
    chkpt_freq = 1.
    pchkpt_freq = 10.

    timestepping_params = TimesteppingParams(
        niter0 = niter0,
        dt = dt,
        ntimesteps_velocity_update = 1000,
        end_time = end_time,
        chkpt_freq = chkpt_freq,
        pchkpt_freq = pchkpt_freq,
        step_thickness = step_thickness,
        verbose = true)

    #output parameters
    # Include all fields needed for forward runs, matching inversion outputs
    outputs = (h = model.fields.gh.h,
            u = model.fields.gu.u,                    # U-grid velocity (for initial_u_veloc)
            uh = model.fields.gh.u,                   # H-grid u velocity (depth-averaged)
            v = model.fields.gv.v,                    # V-grid velocity (for initial_v_veloc)
            vh = model.fields.gh.v,                   # H-grid v velocity (depth-averaged)
            grounded_fraction = model.fields.gh.grounded_fraction,
            beta = model.fields.gh.β,                 # For calculating weertman_c
            preBfactor = model.fields.gh.preBfactor,
            damage = model.fields.g3d.Φ,              # 3D damage field
            eta_av = model.fields.gh.ηav,             # Depth-averaged viscosity
            eta = model.fields.g3d.η,                 # Full 3D viscosity (for initial_viscosity)
            bed_speed = model.fields.gh.bed_speed,    # Needed for weertman_c calculation
            haf = model.fields.gh.haf,                # Height above flotation, needed for weertman_c calculation
            temperature = model.fields.g3d.θ,         # 3D temperature field (for initial_temperature)
            b = model.fields.gh.b,                    # Bed elevation (static, but useful for reference)
            history = model.fields.g3d.strain_history 
            )

    output_freq = 1.
    output_params = OutputParams(outputs = outputs,
                    output_freq = output_freq,
                    output_format = "mat",
                    dump_vel = true,
                    zip_format = "nc",
                    output_path = "./Greenland_CE_600m_relaxation_out/")

    # assemble the simulation
    simulation = Simulation(model = model, 
                        timestepping_params = timestepping_params, 
                        output_params = output_params)

    run_simulation!(simulation)

end

main()