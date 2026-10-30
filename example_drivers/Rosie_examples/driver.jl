#This is a test inversion file to see if WAVI can run  

using WAVI 
using Printf
using ImageFiltering
using SparseArrays

function driver()

#
#Grid and boundary conditions
#
nx = 398
ny = 468
nσ = 26
x0 = -1789000.0
y0 = -833000.0
dx = 2000.0
dy = 2000.0

h_mask=Array{Float64}(undef,nx,ny);
read!("filepath/h_mask.bin",h_mask)
h_mask.=ntoh.(h_mask)

u_iszero=Array{Float64}(undef,nx+1,ny);
read!("filepath/u_iszero.bin",u_iszero)
u_iszero.=ntoh.(u_iszero)

v_iszero=Array{Float64}(undef,nx,ny+1);
read!("filepath/v_iszero.bin",v_iszero)
v_iszero.=ntoh.(v_iszero)

sigma_grid=Array{Float64}(undef,nσ);
read!("filepath/sigma_grid.bin",sigma_grid)
sigma_grid.=ntoh.(sigma_grid)

#Quick fix for sigma
#sigma_grid = ones(size(sigma_grid))-reverse(sigma_grid)

basin_ID=Array{Float64}(undef,nx,ny);
read!("filepath/basinID.bin",basin_ID)
basin_ID.=ntoh.(basin_ID)


grid = Grid(nx = nx, 
            ny = ny,   
            nσ = nσ, 
            x0 = x0, 
            y0 = y0, 
            dx = dx, 
            dy = dy,
            h_mask = h_mask, 
            u_iszero = u_iszero, 
            v_iszero = v_iszero,
            σ = sigma_grid,
            basin_ID = basin_ID)
#
#Bed 
#
bed=Array{Float64}(undef,nx,ny);
read!("filepath/bed.bin",bed)
bed.=ntoh.(bed)

h=Array{Float64}(undef,nx,ny);
read!("filepath/thickness.bin",h)
h.=ntoh.(h) 

 temp=Array{Float64}(undef,nx,ny,nσ);
 read!("filepath/temps.bin",temp)
 temp.=ntoh.(temp)

 #same as accumulation data? Redundant??
 accumulation_data=Array{Float64}(undef,nx,ny);
 read!("filepath/accumulation_data.bin",accumulation_data)
 accumulation_data.=ntoh.(accumulation_data)

#This is from data
 dhdt_data=Array{Float64}(undef,nx,ny);
 read!("filepath/dhdt_data.bin",dhdt_data)
 dhdt_data.=ntoh.(dhdt_data)


initial_conditions = InitialConditions(initial_thickness = h,
#                                        initial_viscosity = viscosity,
                                        initial_temperature = temp)
 #                                       initial_damage = damage)
#                                        initial_u_veloc = gu_u,
 #                                       initial_v_veloc = gv_v)

#we do not evolve the shelves:
evolveShelves = false

accumulation_rate = accumulation_data

params = Params(accumulation_rate = accumulation_rate,
                weertman_c = 0.0,
                evolveShelves = evolveShelves)


 #build the model to use in the initialisation:
println("Starting to make the model")
model = Model(grid = grid,
              bed_elevation = bed, 
              params = params, 
           #   solver_params = solver_params,
              initial_conditions= initial_conditions) 
#
println("The model is made")

println("The inversion test is about to be done")
#perform the inversion

#First read in the data to be used for the inversion:

dhdtaccmask=Array{Float64}(undef,nx,ny);
read!("filepath/dhdt_acc_mask.bin",dhdtaccmask)
dhdtaccmask.=ntoh.(dhdtaccmask)

udata=Array{Float64}(undef,nx+1,ny);
read!("filepath/udata.bin",udata)
udata.=ntoh.(udata)

vdata=Array{Float64}(undef,nx,ny+1);
read!("filepath/vdata.bin",vdata)
vdata.=ntoh.(vdata)

udatamask_raw=Array{Float64}(undef,nx+1,ny);
read!("filepath/udata_mask.bin",udatamask_raw)
udatamask_raw.=ntoh.(udatamask_raw)

vdatamask_raw=Array{Float64}(undef,nx,ny+1);
read!("filepath/vdata_mask.bin",vdatamask_raw)
vdatamask_raw.=ntoh.(vdatamask_raw)

#process the raw data for use in the model:
   udatamask_raw= convert(Array{Bool,2}, udatamask_raw)
   vdatamask_raw= convert(Array{Bool,2}, vdatamask_raw)

   speed_u, udatamask = WAVI.smooth_veloc_data(model.fields.gu, model.fields.gh, udata, udatamask_raw, model.grid.u_isfixed)
   println("speed_u and the udatamask have been returned!")

   speed_v, vdatamask = WAVI.smooth_veloc_data(model.fields.gv, model.fields.gh, vdata, vdatamask_raw, model.grid.v_isfixed)
   println("speed_v and the vdatamask have been returned!")
 

    #check u and v data masks are boolean:
    udatamask_combo = convert(Array{Bool,2}, udatamask)
    vdatamask_combo = convert(Array{Bool,2}, vdatamask)


    dhdtaccmask_combo = ((dhdtaccmask .== 1) .& (model.fields.gh.mask .== 1))
    dhdtaccmask_combo = convert(Array{Bool,2},dhdtaccmask_combo)


    data = DataFields(grid = deepcopy(model.grid),
                      speed_u = speed_u,
                      speed_u_mask =  udatamask_combo,
                      speed_v = speed_v,
                      speed_v_mask =  vdatamask_combo,
                      dhdt = dhdt_data,
                      accumulation_rate = accumulation_rate,
                      dhdtacc_mask=dhdtaccmask_combo)


reltol= 0.5
abstol= 0.1
maxiter= 500
gmres_restart = 50
βgrounded_start=1.e4
βfloating_start=1.0e-10
ηstart_guess = 1.0e7
βpower = 0.2
Bpower_shelf = 0.05
Bpower_grounded = 0.0
#inner_tol = 1.0e-4
#inner_maxiters = 500
cg = false
gmres=true

inversion_params = InversionParams(reltol = reltol,
				    abstol = abstol,
                                    maxiter = maxiter,
                                    gmres_restart = gmres_restart,
                                    βgrounded_start = βgrounded_start,
                                    βfloating_start = βfloating_start,
                                    ηstart_guess = ηstart_guess,
                                    βpower = βpower,
				    Bpower_shelf = Bpower_shelf,
			            Bpower_grounded = Bpower_grounded,
				    cg=cg,
				    gmres=gmres,
				    verbose = true)

#JKVstepping parameters
niter0 = 0
n_iter_out=1.0
max_JKV_iterations = 25
#n_iter_chkpt = 100
n_iter_pchkpt= 5.0

JKVstepping_params = JKVsteppingParams(niter0 = niter0, 
                                     #   n_iter_chkpt = n_iter_chkpt,
                                        n_iter_pchkpt = n_iter_pchkpt,
                                        n_iter_total = max_JKV_iterations,
                                        n_iter_out = n_iter_out)

inversion_output = InversionOutput()
                
   dirichlet_model = deepcopy(model)

   inversion = Inversion(model = dirichlet_model,
                        data = data,
                        inversion_params=inversion_params,
                        inversion_output = inversion_output)

println("About to make inversion_simulation")

outputs = (h = model.fields.gh.h,
            u = model.fields.gu.u,
            uh = model.fields.gh.u,
            v = model.fields.gv.v,
            vh = model.fields.gh.v,
            b = model.fields.gh.b,
            s = model.fields.gh.s,
	    surf_speed=model.fields.gh.surf_speed,
            av_speed = model.fields.gh.av_speed,
            h_mask = model.fields.gh.mask,
            u_mask =  model.fields.gu.mask,
            v_mask =  model.fields.gv.mask,
            grounded_fraction = model.fields.gh.grounded_fraction,
            beta = model.fields.gh.β,
            preBfactor = model.fields.gh.preBfactor,
            JRMS=inversion.inversion_output.JRMS
            )

output_freq = 5
output_params = OutputParams(outputs = outputs, 
                           output_freq = output_freq,
                           output_format = "mat",
                           zip_format = "nc")
   

inversion_simulation = InversionSimulation(model = model, 
                                        inversion = inversion,
                                        JKVstepping_params = JKVstepping_params,
                                        output_params = output_params)
           

#Then run the inversion:
run_inversion_simulation!(inversion_simulation)

#run_inversion!(simulation,inversion)

println("The inversion test has been run")


return model

end