### Zhang Circle of Willis Implementation
# For Modelathon 2026, missing some flow definitions


# ==============================================================================
# Packages
# ==============================================================================

using DifferentialEquations, Plots, Interpolations, Statistics, LinearAlgebra, XLSX, DataFrames, JLD2

# ==============================================================================
# FUNCTIONS
# ==============================================================================

# read parameter values from "ValuesNew.xlsx" for input to model
function ReadValues(filepath)
    # 1. Open the Excel file
    xf = XLSX.readxlsx(filepath)

    # 2. Get 'Sheet1' (default in Excel - change if needed)
    sh = xf["Sheet1"]

    # 3. Set the Data Range manually
    # Set this to read in VALUES ONLY, code not set up to handle header text
    # This returns a generic Matrix of Any type.
    raw_data = sh["M4:P19"]

    # 2. Extract Resistance (R)
    # In the Excel file, this corresponds to Column M
    R = Float64.(raw_data[:, 1])

    # 3. Extract Compliance (C)
    # In Excel: 'C' is Column O, 'C_flag' is Column P
    raw_C = Float64.(raw_data[:, 3])
    c_flags = Float64.(raw_data[:, 4])

    # Apply the flag
    C_values = raw_C .* c_flags

    # 4. Filter out zeros
    C_final = filter(!iszero, C_values)

    return R, C_final
end


# Take moving average
function simple_movmean(data, window_size)
    n = length(data)
    result = zeros(n)
    half_win = floor(Int, window_size / 2)
    for i in 1:n
        # Determine bounds, clamping to edges
        low = max(1, i - half_win)
        high = min(n, i + half_win)
        result[i] = mean(data[low:high])
    end
    return result
end


# interpolate waveforms to provide boundary conditions for model
function define_BCs(t, interpolants)
    # Find relative time within beat (assuming 0.8s period)
    # MATLAB: t = t - 0.8*floor(t/0.8)
    t_rel = mod(t, 0.8)

    Pout = 20.0

    # Use the interpolation objects we created earlier
    PBA = interpolants.PBA_itp(t_rel)
    PICA = interpolants.PICA_itp(t_rel)
    dPBA_dt = interpolants.dPBA_itp(t_rel)

    return PBA, dPBA_dt, PICA, Pout
end


# calculate odes for pressure in model
function cow_zhang_v1_press!(dx, x, p, t)
    # Unpack parameters
    R, C, bc_interpolants = p

    # Map R vector to variables (1-based indexing matches MATLAB)
    R1, R1p, R3, R9, R16 = R[1], R[2], R[3], R[4], R[5]
    R2, R11, R8, R15, R4p = R[6], R[7], R[8], R[9], R[10]
    R4, R6, R13, R7, R5, R12 = R[11], R[12], R[13], R[14], R[15], R[16]

    # Map C vector
    C2, C11, C4p, C4, C7, C5, C12 = C[1], C[2], C[3], C[4], C[5], C[6], C[7]

    # Get Boundary Conditions
    PBA, dPBA_dt, PICA, Pout = define_BCs(t, bc_interpolants)

    # Assign pressures from states
    P2, P4, P5, P7, P11, P12, P14 = x[1], x[2], x[3], x[4], x[5], x[6], x[7]

    # --- Ohm's Law ---
    iBA = PBA / R3
    i4p = (P4 - P5) / R4p
    i4 = (P4 - P12) / R4
    i6 = (P5 - Pout) / R6
    i13 = (P12 - Pout) / R13
    i5 = (P5 - P2) / R5
    i12 = (P12 - P11) / R12
    iICA1 = (PICA - P2) / R1
    iICA2 = (PICA - P11) / R1p
    i9 = (P2 - Pout) / R9
    i16 = (P11 - Pout) / R16
    i2 = (P2 - P7) / R2
    i11 = (P11 - P14) / R11
    i8 = (P7 - Pout) / R8
    i7 = (P7 - P14) / R7
    i15 = (P14 - Pout) / R15

    # --- Continuity ---
    # iC4p = i4p - (i06+ i5)

    # --- Derivatives ---
    dx[3] = (1 / R4p * P4 + 1 / R5 * P2 - (1 / R4p + 1 / R6 + 1 / R5) * P5) / C4p # dP5

    dx[4] = (1 / R2 * P2 + 1 / R7 * P14 - (1 / R2 + 1 / R8 + 1 / R7) * P7) / C2   # dP7

    dx[6] = (1 / R4 * P4 + 1 / R12 * P11 - (1 / R4 + 1 / R13 + 1 / R12) * P12) / C4 # dP12

    dx[1] = (1 / R5 * P5 + 1 / R2 * P7 + 1 / R1 * PICA - (1 / R1 + 1 / R5 + 1 / R9 + 1 / R2) * P2) / C5 # dP2

    dx[5] = (1 / R11 * P14 + 1 / R1p * PICA + 1 / R12 * P12 - (1 / R12 + 1 / R15 + 1 / R11 + 1 / R1p) * P11) / C12 # dP11

    dx[7] = (1 / R7 * P7 + 1 / R11 * P11 - (1 / R7 + 1 / R15 + 1 / R11) * P14) / (C11 + C7) # dP14

    # Complex expression for dP4
    num_dP4 = R4 * R4p * dPBA_dt + R3 * R4p * dx[6] + R3 * R4 * (iC4p / C4p)
    den_dP4 = R4 * R4p + R3 * R4p + R3 * R4
    dx[2] = num_dP4 / den_dP4 # dP4
end


# calculate flows from rpessure values
function cow_zhang_v1_flows(t_vec, x_mat, p)
    R, C, bc_interpolants = p

    # Map Parameters
    R1, R1p, R3, R9, R16 = R[1], R[2], R[3], R[4], R[5]
    R2, R11, R8, R15, R4p = R[6], R[7], R[8], R[9], R[10]
    R4, R6, R13, R7, R5, R12 = R[11], R[12], R[13], R[14], R[15], R[16]
    C2, C11, C4p, C4, C7, C5, C12 = C[1], C[2], C[3], C[4], C[5], C[6], C[7]

    # Extract States
    P2, P4, P5, P7, P11, P12, P14 = [x_mat[:, i] for i in 1:7]

    # Calculate BCs for all time points
    bcs = define_BCs.(t_vec, Ref(bc_interpolants))
    PBA = [b[1] for b in bcs]
    PICA = [b[3] for b in bcs]
    Pout = 20.0 # From define_BCs

    ## Ohm's law
    iBA = PBA ./ R3
    i4p = (P4 .- P5) ./ R4p 
    i4 = (P4 .- P12) ./ R4
    i6 = (P5 .- Pout) ./ R6
    i13 = (P12 .- Pout) ./ R13
    i5 = (P5 .- P2) ./ R5
    i12 = (P12 .- P11) / R12
    iICA1 = (PICA .- P2) ./ R1
    iICA2 = (PICA .- P11) ./ R1p
    i9 = (P2 .- Pout) ./ R9
    i16 = (P11 .- Pout) / R16
    i2 = (P2 .- P7) ./ R2
    i11 = (P11 .- P14) ./ R11
    i8 = (P7 .- Pout) ./ R8
    i7 = (P7 .- P14) ./ R7
    i15 = (P14 .- Pout) ./ R15
    ## Continuity equations 
    i4pp = i6 .+ i5
    iC4p = i4p .- i4pp
    i2p = i8 .+ i7
    iC2 = i2 .- i2p
    i4z = i12 .+ i13
    iC4 = i4 .- i4z
    i5p = i9 .+ i2 .- iICA1
    iC5 = i5 .- i5p
    i12p = i16 .+ i11 .- iICA2
    iC12 = i12 .- i12p


    # Derivatives for missing flows
    dt_vec = diff(t_vec)
    dP14_diff = diff(P14)

    # Calculate iC7
    iC7_raw = (dP14_diff ./ dt_vec) .* C7
    iC7 = [iC7_raw; 0.0]

    # Calculate iC11
    iC11_raw = (dP14_diff ./ dt_vec) .* C11
    iC11 = [iC11_raw; 0.0]

    i11p = i11 .- iC11
    i7p = i7 .- iC7

    # Concatenate into Q matrix (Columns match MATLAB order)
    q = hcat(iBA, i4p, i4pp, iC4p, i6, i5, i5p, iC5, iICA1, i9, i2, i2p, iC2, i8, i7, i7p, iC7, i4, i4z, iC4, i13, i12, i12p, iC12, iICA2, i16, i11, i11p, iC11, i15)
    return q
end



# ==============================================================================
# Main Code
# ==============================================================================


## 1. PRE-PROCESSING (Signal Interpolation)
# ==============================================================================

# Define the names of Excel files containing values and input waveforms
valfile = joinpath(@__DIR__, "ValuesNew.xlsx")

# Use functions to read in params and waveforms
R, C = ReadValues(valfile)

# Load waveforms
temp = load(joinpath(@__DIR__, "Zhang_digitise.jld2"))
t_ref_raw, BAx, BAy, ICAx, ICAy = temp["waveforms"]

# Interpolate pressures onto chosen time vector (t_ref_raw)
itp_BA = LinearInterpolation(BAx, BAy, extrapolation_bc=Line())
itp_ICA = LinearInterpolation(ICAx, ICAy, extrapolation_bc=Line())

BAy_ref = itp_BA.(t_ref_raw)
ICAy_ref = itp_ICA.(t_ref_raw)

# Apply periodic moving average for C2 continuity
N_smooth = 7
k = floor(Int, N_smooth / 2)

# Wrap array for periodic smoothing: [End bit; Middle; Start bit]
wk1_BA = [BAy_ref[end-k+1:end]; BAy_ref; BAy_ref[1:k]]
wk1wrap_BA = simple_movmean(wk1_BA, N_smooth)
BAy_ref = wk1wrap_BA[k+1:end-k]

wk1_ICA = [ICAy_ref[end-k+1:end]; ICAy_ref; ICAy_ref[1:k]]
wk1wrap_ICA = simple_movmean(wk1_ICA, N_smooth)
ICAy_ref = wk1wrap_ICA[k+1:end-k]

# Compute forward difference for PBA
wk1_diff = diff(BAy_ref)
# Append end difference to maintain length
push!(wk1_diff, (BAy_ref[1] - BAy_ref[end]))
wk1_diff = wk1_diff ./ mean(diff(t_ref_raw))
dPBA_dt_ref = wk1_diff

# Create Fast Interpolation Objects for the ODE solver
# We wrap these in a NamedTuple to pass to the solver efficiently
bc_interpolants = (
    PBA_itp=LinearInterpolation(t_ref_raw, BAy_ref, extrapolation_bc=Periodic()),
    PICA_itp=LinearInterpolation(t_ref_raw, ICAy_ref, extrapolation_bc=Periodic()),
    dPBA_itp=LinearInterpolation(t_ref_raw, dPBA_dt_ref, extrapolation_bc=Periodic())
)

# Plot Pre-processing Check
p_check = Plots.plot(layout=(1, 2), size=(1000, 400))
Plots.plot!(p_check[1], t_ref_raw, BAy_ref, label="Pressure", ylabel="Pressure")
Plots.plot!(p_check[1], t_ref_raw, dPBA_dt_ref, label="dP/dt", ylabel="dP/dt", secondary=true)
Plots.plot!(p_check[2], t_ref_raw, ICAy_ref, label="ICA", color=:red)
display(p_check)




## 2. SOLVER
# ==============================================================================

# Simulation Extent
t0 = 0.0
tf = 9 * 0.8
tspan = (t0, tf)

# Initial Values
mcfp = 40.0
x0 = fill(mcfp, 7)

# Pack Parameters
p = (R, C, bc_interpolants)

# Solver
# Note: ode15s (MATLAB) is for stiff systems. TRBDF2() or Rosenbrock23() are good Julia equivalents.
prob = ODEProblem(cow_zhang_v1_press!, x0, tspan, p)
println("Solving System...")
@time sol = solve(prob, TRBDF2(), dtmax=0.01)

# Extract Data
t_arr = sol.t
x_mat = stack(sol.u)'

# Calculate Flows
println("Calculating Flows...")
q = cow_zhang_v1_flows(t_arr, x_mat, p)



## 3. PLOT RESULTS
# ==============================================================================

press_legend1 = ["ICA L out (P2)" "BA out (P4)" "PCA2 L in (P5)" "ACA2 L in (P7)" "ICA R out (P11)" "PCA2 R in (P12)" "ACA2 R out (P14)"] # 1xN Matrix
press_legend3 = ["02", "04", "05", "07", "11", "12", "14"] # Vector for loop
flow_legend = ["BA (iBA)" "PCA1 L in (i4p)" "PCA1 L out (i4pp)" "PCA1 L stored (iC4p)" "PCA2 L (i6)" "PCoA L in (i5)" "PCoA L out (i5p)" "PCoA L stored (iC5)" "ICA L (iICA1)" "MCA L (i9)" "ACA1 L in (i2)" "ACA1 L out (i2p)" "ACA1 L stored (iC2)" "ACA2 L (i8)" "ACoA in (i7)" "ACoA out (i7p)" "ACoA stored (iC7)" "PCA1 R in (i4)" "PCA1 R out (i4z)" "PCA1 R stored (iC4)" "PCA2 R (i13)" "PCoA R in (i12)" "PCoA R out (i12p)" "PCoA R stored (iC12)" "ICA R (iICA2)" "MCA R (i16)" "ACA1 R in (i11)" "ACA1 R out (i11p)" "ACA1 R stored (iC11)" "ACA2 R (i15)"]

# Select cycle (Last 2 cycles approx)
N_cycles = 1
indices = findall(t -> t >= t_arr[end] - N_cycles * 0.8, t_arr)

# 1. Summary Pressure Plot
p_summary = Plots.plot(t_arr, x_mat, label=press_legend1, xlabel="Time (s)", ylabel="Pressure (mmHg)", title="Pressures Summary")
display(p_summary)

# 2. Separated Pressures
# We use a layout to mimic the subplot loop
p_separated = Plots.plot(layout=(7, 1), size=(600, 1000), link=:x)

for i in 1:7
    # Note: Reverse order 8-i to match MATLAB's top-down fill
    Plots.plot!(p_separated[8-i], t_arr[indices], x_mat[indices, i],
        label=press_legend1[i],
        ylabel="P (mmHg)",
        ylim=(40, 130))
end
Plots.xlabel!(p_separated[7], "Time (s)")
display(p_separated)

# 3. Flows (Optional - uncomment to view)
# flow_legend = ["iBA" "i4p" ... ] # Define all 30 strings here if needed

p_flows = Plots.plot(layout=(5, 6), size=(1200, 800), link=:x)


for i in 1:30
    Plots.plot!(p_flows[i], t_arr[indices], q[indices, i],
        label=flow_legend[i],
        ylabel="Q (ml.s-1)",
        ylim=(-2, 10))
end
Plots.xlabel!(p_flows[30], "Time (s)")
display(p_flows)


display(p_flows)


# Save the summary plot
savefig(p_summary, "pressures_summary.png")

# Save the separated subplot layout
savefig(p_separated, "pressures_separated.png")

# Add dpi=300 for high-quality (publication ready) output
p_summary = Plots.plot(t_arr, x_mat, dpi=300)
savefig(p_summary, "high_quality_summary.png")

