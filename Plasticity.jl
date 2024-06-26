using Waveforms
using Plots

# Helper function to get a triangle wave strain time history
function get_cycles(; peak2peak=1, num_cycles=1, num_steps_per_cycle=100)
    return peak2peak * trianglewave1.(range(0.0, stop=num_cycles, length=num_steps_per_cycle*num_cycles))
end

ϵ = get_cycles(peak2peak=0.3/100, num_cycles=1, num_steps_per_cycle=1000)
σₒ = 250.0
E = 200000.0
H_iso = 10000.0
H_kin = 40000.0

# Function for isotropic hardening
function Isotropic_Hardening(ϵ, σₒ, E, H_iso)
    σ = zeros(length(ϵ))
    εₚ = 0.0
    σᵧ = σₒ # Yield stress

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

function Kinematic_Hardening(ϵ, σₒ, E, H_kin)
    σ = zeros(length(ϵ))
    εₚ = 0.0
    σᵧ = σₒ # Yield stress
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

function Combined_Hardening(ϵ, σₒ, E, H_iso, H_kin)
    σ = zeros(length(ϵ))
    εₚ = 0.0
    σᵧ = σₒ # Yield stress
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

σ_iso = Isotropic_Hardening(ϵ, σₒ, E, H_iso)
σ_kin = Kinematic_Hardening(ϵ, σₒ, E, H_kin)
σ_com = Combined_Hardening(ϵ, σₒ, E, H_iso, H_kin)
σ_com_pla = Combined_Hardening(ϵ, σₒ, E, 0.0, 0.0)
plot(ϵ, σ_com_pla, label="Perfect Plasticity", xlabel="ϵ", ylabel="σ", title="Plot of σ vs ϵ")
plot!(ϵ, σ_iso, label="Isotropic", xlabel="ϵ", ylabel="σ", title="Plot of σ vs ϵ")
plot!(ϵ, σ_kin, label="Kinematic", xlabel="ϵ", ylabel="σ", title="Plot of σ vs ϵ")
plot!(ϵ, σ_com, label="Combined", xlabel="ϵ", ylabel="σ", title="Plot of σ vs ϵ")