# module Analysis
#
# using PyCall, Plots#; pyplot()
#
# ase_read = pyimport("ase.io")["read"]
# ase_write = pyimport("ase.io")["write"]
#
# export plotV, ase_read, ase_write
#
# function plotV(al :: Array{PyCall.PyObject,1}; hist=true)
#     Vl = []
#
#     for at in al
#         V = dot(at[:cell][3,:],
#         cross(at[:cell][1,:], at[:cell][2,:]))
#         push!(Vl, V)
#     end
#
#     if hist == true
#         histogram(Vl)
#     else
#         plot(Vl)
#     end
# end
#
# end
