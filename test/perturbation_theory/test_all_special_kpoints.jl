# Comprehensive test: frequency_shifts for every special k-point across all plane groups
# and space groups, with and without time-reversal (realify), both polarizations (2D),
# and degeneracy_idx ∈ {1, 2}.
#
# This test is intentionally broad: it verifies that the perturbation theory machinery
# can construct a model without erroring. Failures are logged to a text file rather than
# causing the test suite to abort, so we can track remaining issues systematically.

using Test, Crystalline, LinearAlgebra
using EmptyLattice
using EmptyLattice.PerturbationTheory

# Log file for failures — written next to runtests.jl
enable_logging = false # Set to true to write a log file of all failures
LOGFILE = joinpath(@__DIR__, "..", "frequency_shifts_failures.log")

function _run_all_special_kpoints()
    failures = []  # (D, sgnum, klab, polarization, timereversal, idx, err)
    successes = 0
    for D in (2, 3)
        polarizations = D == 2 ? (:TM, :TE) : (nothing,)
        for sgnum in 1:MAX_SGNUM[D]
            cntr = centering(sgnum, D)
            Rs = primitivize(directbasis(sgnum, Val(D)), cntr)
            Gs = dualbasis(Rs)
            for timereversal in (false, true)
                # Get all irreps; optionally realify for time-reversal
                lgirsd = lgirreps(sgnum, Val(D))
                timereversal && realify!(lgirsd)
                lgirsd = primitivize(lgirsd)
                for (klab, lgirs) in lgirsd
                    # Skip non-special k-points (have free parameters)
                    isspecial(first(lgirs)) || continue
                    for polarization in polarizations
                        for idx in (1, 2)
                            idx == 1 && klabel(lgirs) == "Γ" && continue  # skip q = 0 case
                            try
                                frequency_shifts(lgirs, Gs, idx; polarization)
                                successes += 1
                                msg = ""
                            catch e
                                msg = sprint(showerror, e)
                                # Truncate long error messages
                                if length(msg) > 200
                                    msg = msg[1:200] * "..."
                                end
                                push!(failures, (D, sgnum, klab, polarization, timereversal, idx, msg))
                            end
                        end
                    end
                end
            end
        end
    end

    return successes, failures
end

@testset "All special k-points: frequency_shifts smoke test" begin
    successes, failures = _run_all_special_kpoints()

    # Write failures to log file
    if enable_logging
        open(LOGFILE, "w") do io
            println(io, "# frequency_shifts failure log")
            println(io, "# Generated: ", round(Int, time()))
            println(io, "# Format: D | sgnum | klab | polarization | time_reversal | degeneracy_idx | error")
            println(io, "#")
            if isempty(failures)
                println(io, "# No failures!")
            else
                println(io, "# Total failures: ", length(failures))
                println(io, "# Total successes: ", successes)
                println(io)
                for (D, sgnum, klab, pol, tr, idx, msg) in failures
                    tr_str = tr ? "TR" : "noTR"
                    println(io, "$D | $sgnum | $klab | $pol | $tr_str | $idx | $msg")
                end
            end
        end
    end

    # Report
    n_fail = length(failures)
    @info "frequency_shifts smoke test: $successes succeeded, $n_fail failed"
    if n_fail > 0
        enable_logging && @info "Failure log written to $LOGFILE"
        # Show a summary of unique error types
        err_types = Dict{String, Int}()
        for (_, _, _, _, _, _, msg) in failures
            # Extract first line of error as a key
            key = first(split(msg, '\n'))
            err_types[key] = get(err_types, key, 0) + 1
        end
        @info "Error type summary:" err_types
    end

    # We expect successes > 0; the test records failures rather than failing hard
    @test successes > 0
    # Record the failure count so we can track regressions / improvements over time.
    # Currently fails at P/PA points of I-centered cubic groups (199, 206, 214, 220, 230)
    # This is essentially Crystalline.jl's issue #12, so it cannot be directly fixed here.
    @test n_fail == 14  
end
