using Waveforms
using Plots

function get_cycles(; peak2peak=1, num_cycles=1, num_steps_per_cycle=100)
    return peak2peak * trianglewave1.(range(0.0, stop=num_cycles, length=num_steps_per_cycle*num_cycles))
end

ϵ = get_cycles(peak2peak=0.3/100, num_cycles=1, num_steps_per_cycle=1000)
σₒ = 250.0
E = 200.0e3
H_iso = 10.0e3
H_kin = 40.0e3

# Function for isotropic hardening
function Isotropic_Hardening(ϵ, σₒ, E, H_iso)
    σ = zeros(length(ϵ))
    εₚ = 0.0
    σᵧ = σₒ 

    for i in 1:length(ϵ)
        εₑ = ϵ[i] - εₚ
        σ_trial = E * εₑ
        f = abs(σ_trial) - σᵧ
        if f <= 0
            σ[i] = σ_trial
        else
            Δλ = f / (E + H_iso)
            εₚ += Δλ * sign(σ_trial)
            σᵧ += H_iso * Δλ
            σ[i] = E * (ϵ[i] - εₚ)
        end
    end
return σ
end

# Function for kinematic hardening
function Kinematic_Hardening(ϵ, σₒ, E, H_kin)
    σ = zeros(length(ϵ))
    εₚ = 0.0
    σᵧ = σₒ 
    α = 0.0

    for i in 1:length(ϵ)
        εₑ = ϵ[i] - εₚ
        σ_trial = E * εₑ
        f = abs(σ_trial - α) - σᵧ
        if f <= 0
            σ[i] = σ_trial
        else
            Δλ = f / (E + H_kin)
            εₚ += Δλ * sign(σ_trial - α)
            α += H_kin * Δλ * sign(σ_trial - α)
            σ[i] = E * (ϵ[i] - εₚ)
        end
    end
return σ
end

# Function for combined isotropic and kinematic hardening
function Combined_Hardening(ϵ, σₒ, E, H_iso, H_kin)
    σ = zeros(length(ϵ))
    εₚ = 0.0
    σᵧ = σₒ 
    α = 0.0

    for i in 1:length(ϵ)
        εₑ = ϵ[i] - εₚ
        σ_trial = E * εₑ
        f = abs(σ_trial - α) - σᵧ
        if f <= 0
            σ[i] = σ_trial
        else
            Δλ = f / (E + H_kin + H_iso)
            εₚ += Δλ * sign(σ_trial - α)
            α += H_kin * Δλ * sign(σ_trial - α)
            σ[i] = E * (ϵ[i] - εₚ)
            σᵧ += H_iso * Δλ
        end
    end
return σ
end

σ_pp = Combined_Hardening(ϵ, σₒ, E, 0.0, 0.0)
σ_iso = Combined_Hardening(ϵ, σₒ, E, H_iso, 0.0)
σ_kin = Combined_Hardening(ϵ, σₒ, E, 0.0, H_kin)
σ_com = Combined_Hardening(ϵ, σₒ, E, H_iso, H_kin)
plot(ϵ, σ_pp, label="Perfect Plasticity", xlabel="ϵ", ylabel="σ", title="σ vs ϵ")
plot!(ϵ, σ_iso, label="Isotropic", xlabel="ϵ", ylabel="σ", title="σ vs ϵ")
plot!(ϵ, σ_kin, label="Kinematic", xlabel="ϵ", ylabel="σ", title="σ vs ϵ")
plot!(ϵ, σ_com, label="Combined", xlabel="ϵ", ylabel="σ", title="σ vs ϵ")