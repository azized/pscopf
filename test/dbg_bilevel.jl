using JuMP, BilevelJuMP, Cbc

#=
    TS: [11h, 11h15]
    S: [S1,S2]
                        bus 1                   bus 2
                        |                      |
    (limitable) wind_1_1|       "1_2"          |
    Pmin=0, Pmax=PMAX1  |                      |
    Csta=0, Cprop=1     |                      |
      P1                |----------------------|
                        |        L12           |
                        |                      |
                        |                      |
           load(bus_1)  |                      |load(bus_2)
      D1                |                      | D2
=#


model = BilevelModel(Cbc.Optimizer)
#=
BilevelJuMP.AbstractBoundedMode{T}
BilevelJuMP.ComplementMode{T}
BilevelJuMP.IndicatorMode{T}
BilevelJuMP.NoMode{T}
BilevelJuMP.ProductMode{T}
BilevelJuMP.SOS1Mode{T}
BilevelJuMP.StrongDualityMode{T}
=#
BilevelJuMP.set_mode(model, BilevelJuMP.IndicatorMode())

PMAX1 = 100.

# equivalent testcase : no_problem
# works as expected
P1 = 40.
D1 = 10.
D2 = 30.
L12 = 35.

# # equivalent testcase : EOD_problem_needs_capping
# # works as expected
# P1 = 100.
# D1 = 10.
# D2 = 30.
# L12 = 35.

# # equivalent testcase : EOD_problem_needs_loss_of_load
# # FAIL : cuts all prod and conso for a cost of 130*e17+40
# # while a solution of cost 90*1e4 exists
# # adding @constraint(Lower(model), e_st == 0), gives the expected result
# P1 = 40.
# D1 = 100.
# D2 = 30.
# L12 = 35.

# # equivalent testcase : RSO_problem
# # FAIL : cuts all prod and conso for a cost of 50*e17+50
# # while a solution of cost 5*1e4+5 exists
# # adding @constraint(Lower(model), e_st == 5), gives the expected result
# P1 = 50.
# D1 = 10.
# D2 = 40.
# L12 = 35.


#upper
@variable(Upper(model), lolmin_st, lower_bound=0., upper_bound=D1+D2)
@variable(Upper(model), lol_n1st, lower_bound=0., upper_bound=D1)
@variable(Upper(model), lol_n2st, lower_bound=0., upper_bound=D2)

@variable(Upper(model), emin_st, lower_bound=0., upper_bound=P1)
@variable(Upper(model), e_n1st, lower_bound=0., upper_bound=P1)
@variable(Upper(model), penr_n1st, lower_bound=0., upper_bound=P1)#modif
@variable(Upper(model), plim_n1st, lower_bound=0., upper_bound=PMAX1)#modif
@variable(Upper(model), islimited_n1st, binary=true)
@variable(Upper(model), plimxislim_n1st, lower_bound=0., upper_bound=PMAX1)#modif

#lower
@variable(Lower(model), lol_st, lower_bound=0., upper_bound=D1+D2)

@variable(Lower(model), e_st, lower_bound=0., upper_bound=P1)
@variable(Lower(model), penr_st, lower_bound=0., upper_bound=P1)


# Upper
@objective(Upper(model), Min, 1e4*lolmin_st+emin_st)
@constraints(Upper(model), begin
    -35. <= 0.5*penr_n1st +0.5*lol_n1st -0.5*D1 -0.5*lol_n2st +0.5*D2 <= 35. #RSO
    lol_st == lol_n2st + lol_n1st
    e_n1st == P1 - penr_n1st
    #plimxislim_n1st = plim_n1st * islimited_n1st
    plimxislim_n1st <= plim_n1st
    plimxislim_n1st <= PMAX1*islimited_n1st
    plim_n1st <= PMAX1*(1-islimited_n1st) + plimxislim_n1st
    # penr_n1st min(P1, plim_n1st)
    #penr_n1st <= P1 #as bound
    penr_n1st <= plim_n1st
    penr_n1st == plimxislim_n1st + P1*(1 - islimited_n1st)
    penr_st == penr_n1st
    e_st == e_n1st
end)

# Lower
@objective(Lower(model), Min, 1e4*lol_st+e_st)
@constraints(Lower(model), begin
    P1-e_st == D1 + D2 - lol_st
    emin_st <= e_st
    lolmin_st <= lol_st
end)

optimize!(model)
# optimize!(model, lower_prob="low", upper_prob="up", bilevel_prob="bil", solver_prob="solv")

println("obj lower : ", objective_value(Lower(model)))
println("obj upper : ", objective_value(model))

println("lolmin : ", value(lolmin_st))
println("emin : ", value(emin_st))
println("lol : ", value(lol_st))
println("e : ", value(e_st))
