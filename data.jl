using CSV
using Query
using DataFrames

casesf = reformat(CSV.read("time_series_covid19_confirmed_global.csv"), :Date, :Cases)
deathf = reformat(CSV.read("time_series_covid19_deaths_global.csv"), :Date, :Deaths)
recovf = reformat(CSV.read("time_series_covid19_recovered_global.csv"), :Date, :Recoveries)


chinaf = joindates(select_country.([casesf, deathf, recovf], "China")...)
hubeif = joindates(select_province.([casesf, deathf, recovf], "China", "Hubei")...)

parsedates!(hubeif)

# cdf = @from i in casesf begin
#     @join j in deathf on i.Date equals j.Date #and i.Province equals j.Province
#     @where i.Province == j.Province  && i.Country == j.Country
#     @select {i.Date,i.Cases,j.Deaths,j.Province, j.Country}
#     @collect DataFrame
# end

caseplot(hubeif)

df = @from i in joindates(casesf, deathf, recovf) begin
  @where i.Cases != 0
  @select i
  @collect DataFrame
end

CSV.write("covid.csv", df)


italyf = joindates(select_country.([casesf, deathf, recovf], "Italy")...)
francef = joindates(select_mainland.([casesf, deathf, recovf], "France")...)
ukf = joindates(select_mainland.([casesf, deathf, recovf], "United Kingdom")...)
euf = vcat(italyf, francef, ukf)

CSV.write("eucovid.csv", euf)
