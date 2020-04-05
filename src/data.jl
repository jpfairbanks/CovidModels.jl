module Data
using DataFrames
using Query
using Dates
using StatsPlots

export reformat, select_country, select_mainland, select_province, joindates, parsedates!, caseplot

function reformat(df::DataFrame, variable_name=:Date, value_name=:Count, cols=5)
    df = stack(df,names(df)[cols:end],
               variable_name=variable_name,
               value_name=value_name)
end

"""    select_province(df, country, province)

selects the rows of a dataframe based on their country and province
"""
function select_province(df, country, province)
    df = @from i in df begin
        @where i.Country == country && i.Province == province
        @select i
        @collect DataFrame
    end
end

"""    select_country(df, country, province)

selects the rows of a dataframe based on their country
"""
function select_country(df, country)
    df = @from i in df begin
        @where i.Country == country
        @select i
        @collect DataFrame
    end
end

"""    select_mainland(df, country, province)

selects the rows of a dataframe based on their country and excludes overseas territories
"""
function select_mainland(df, country)
    df = @from i in df begin
        @where i.Country == country && isna(i.Province)
        @select i
        @collect DataFrame
    end
end

"""     joindates(casesf, deaths, recovf)

join the cases, deaths, and recovery frames into one dataframe based on the dates and country+province
"""
function joindates(casesf, deathf, recovf)
    cdf = casesf |> @join(deathf, {_.Date, _.Province, _.Country},
    {_.Date, _.Province, _.Country},
    {_.Date, _.Province, _.Country, _.Cases, __.Deaths}) |> DataFrame

    ff = cdf |> @join(recovf, {_.Date, _.Province, _.Country},
    {_.Date, _.Province, _.Country},
    {_.Date, _.Province, _.Country, _.Cases, _.Deaths, __.Recoveries}) |> DataFrame
    return ff
end

"""     parsedates!(df::DataFrame)

replace the Dates column of a dataframe as string or symbol with the parsed Dates.Date version
"""
function parsedates!(df::DataFrame)
    df.Date = Date.(string.(df.Date), dateformat"m/d/y")
    return df
end

"""    caseplot(df)

plot the number of cases, deaths, recoveries over time
"""
function caseplot(df)
    @df df plot(:Date, [:Cases, :Deaths, :Recoveries], label=["cases" "deaths" "recoveries"])
end
end
