using Colors
using DataFrames
using CSV

# RGB.(colormatch.(350:10:750))

df = DataFrame(CSV.File("catalog.dat",delim=" ",ignorerepeated=true,header=0))