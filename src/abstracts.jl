abstract type AbstractDataGenerator end
function launch(launchable::AbstractDataGenerator) error("unimplemented") end

abstract type  AbstractContext end
abstract type  AbstractSchedule end


"""
    AbstractRunnable

an abstract type interface for all decisional and utilitary processes (TSO, Market and assessment steps).
Such structure might need to overload the following functions :
- `compute_firmness(runnable::AbstractRunnable, ech, next_ech, TS, context::AbstractContext)`
- `run(runnable::AbstractRunnable, ech, firmness, TS, context::AbstractContext)`
- `affects_market_schedule(runnable::AbstractRunnable)`
- `update_market_schedule!(context::AbstractContext, ech, result, firmness, runnable::AbstractRunnable)`
- `affects_tso_schedule(runnable::AbstractRunnable)`
- `update_tso_schedule!(context::AbstractContext, ech, result, firmness, runnable::AbstractRunnable)`
- `affects_tso_actions(runnable::AbstractRunnable)`
- `update_tso_actions!(context::AbstractContext, ech, result, firmness, runnable::AbstractRunnable)`
"""
abstract type  AbstractRunnable end
abstract type  AbstractRunnableConfigs end
"""
Determines the firmnes of commitment and powerlevel decisions for each timestep
"""
function compute_firmness(runnable::AbstractRunnable, ech, next_ech, TS, context::AbstractContext) error("unimplemented") end
"""
Launches the runnable :
    get relevant data from context
    build the optimization model
    solve the model
    return a result object to be exploited by the updating methods
"""
function run(runnable::AbstractRunnable, ech, firmness, TS, context::AbstractContext) error("unimplemented") end
"""
If returns True, will launch the market_schedule related updates during sequence execution
"""
function affects_market_schedule(runnable::AbstractRunnable) return false end
"""
Updates the market schedule using the result of run()
"""
function update_market_schedule!(context::AbstractContext, ech, result, firmness, runnable::AbstractRunnable) error("unimplemented") end
"""
If returns True, will launch the tso_schedule related updates during sequence execution
"""
function affects_tso_schedule(runnable::AbstractRunnable) return false end
"""
Updates the TSO schedule using the result of run()
"""
function update_tso_schedule!(context::AbstractContext, ech, result, firmness, runnable::AbstractRunnable) error("unimplemented") end
"""
If returns True, will launch the tso actions related updates during sequence execution
"""
function affects_tso_actions(runnable::AbstractRunnable) return false end
"""
Updates the TSO actions using the result of run()
"""
function update_tso_actions!(context::AbstractContext, ech, result, firmness, runnable::AbstractRunnable) end
#FIXME are these compatible with the Assessment step ?

abstract type  AbstractTSO <: AbstractRunnable  end
function affects_tso_schedule(runnable::AbstractTSO) return true end
function affects_tso_actions(runnable::AbstractTSO) return true end

abstract type  AbstractMarket <: AbstractRunnable  end
function affects_market_schedule(runnable::AbstractMarket) return true end

abstract type  AbstractAssessment <: AbstractRunnable end
# needs to take decisions => updates a planning/schedule

abstract type DeciderType end
