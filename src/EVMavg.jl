"""
Package EVMavg: Average values by the Expected Value Method (EVM)

### Author(s):
Mike Kretlow <mike@kretlow.de> [github.com/mkretlow]

### Todo:
* EVM: handle symmetric and asymmetric errors

"""
module EVMavg

export evm, avg


"""
    uavg, wavg, se = avg(data, sigmas)

Computes the unweighted and the weighted average and the standard error of the result (wavg). The weights are defined by w = 1/sigma^2

### Arguments
* `data::AbstractArray`: data
* `sigmas::AbstractArray`: sigmas (errors) of the data

### Returns
* Returns the unweighted average `uavg`, the weighted average `wavg` and the stander error `se` of the average

"""
function avg(data::AbstractArray,sig::AbstractArray)

    n = length(data)

    if n == 1 return data[1],data[1],sig[1] end

    w = 1 ./sig.^2
    sw = sum(w)

    avgu = sum(data) / n # unweighted avg
    avgw = 0

    for i = 1:n
        avgw = avgw + data[i]*w[i]
    end

    avgw = avgw / sw    # weighted avg

    ss = 0
    for i = 1:n
        ss = ss + (data[i]-avgu)^2
    end

    ss = ss/(n-1)
    b  = sw^2 / sum(w.^2)
    se = sqrt(ss/b)     # standard error of avgw

    return avgu, avgw, se

end


"""
    evm_avg, se_int, se_ext = evm(data, sigmas)

Computes the EVM average.

### Arguments
* `data::AbstractArray`: data
* `sigmas::AbstractArray`: sigmas (errors) of the data

### Returns
* Returns the EVM average `evm_avg`, and the internal `se_int` and external `se_ext` standard error of the average

### References
This method is described in Birch, M., Singh, B., 2014, Nuclear Data Sheets, 120, 106.

"""
function evm(data::AbstractArray,err::AbstractArray)

    n = length(data)

    n == 1 && return data[1],err[1],err[1]

    M = zeros(n)
    w = zeros(n)

    for i = 1:n
        M[i] = 0
        for j = 1:n
            M[i] = M[i] + evm_n(data[i],data[j],err[j])
            #M[i] = M[i] + evm_a(data[i],data[j],err[j],err[j])
        end
        M[i] = M[i]/n
    end

    for i = 1:n
        w[i] = M[i] / sum(M)
    end

    xb = 0
    for i = 1:n
        xb = xb + w[i]*data[i]
    end

    se_int = 0
    for i = 1:n
        se_int = se_int + w[i]^2*err[i]^2
    end

    se_int = sqrt(se_int)

    se_ext = 0
    for i = 1:n
        se_ext = se_ext + w[i]*(data[i]-xb)^2
    end

    se_ext = sqrt(se_ext)

    return xb, se_int, se_ext

end


# Helper function for asymmetric standard errors
function evm_a(x,μ,a,b)

    if x <= my  res = sqrt(2/(pi*(a+b)^2)) * exp( -(x-μ)^2 / (2*b^2)) end
    if x >  my  res = sqrt(2/(pi*(a+b)^2)) * exp( -(x-μ)^2 / (2*a^2)) end

end


# Helper function for standard symmetric errors
function evm_n(x,μ,ss)

    res = 1/(sqrt(2*pi)*ss) * exp( -(x-μ)^2/(2*ss^2))

end


end # module
