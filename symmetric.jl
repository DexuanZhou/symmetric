using Combinatorics
using BenchmarkTools
using Polynomials
using Plots
using Jacobi
using SparseArrays

##################
#sum_Tk 和 Prod_Tk是为了求laplace
#sum_Tk为(1)(2)
#Prod_Tk为(1)(2)中每一个乘积项
##################
function basis_function(N,d,x,k,range)
    Prod_Tk=zeros(N,N)
    Sum_Tk=zeros(1,N)
        for i=1:N
            for z=1:N
                Prod_temp=1
                for j=1:d
                    Prod_temp=Prod_temp*chebyshev(x[z,j],k[i,j])
                end
                Prod_Tk[i,z]=Prod_temp
            end
            Sum_Tk[i]=sum(Prod_Tk[i,:])
        end
    boundary_func=vcat(map(x->range^2-x^2, x))
    func_value=prod(Sum_Tk)*prod(boundary_func)
    #Sum_Tk=zeros(Float64,N)
    #Prod_Tk=zeros(Float64,N,N)
    #temp=[1.0*chebyshev(x[z,j],k[i,j]) for i=1:N for z=1:N for j=1:d]
    #temp=reshape(temp,d,N*N)
    #Prod_Tk=prod(temp,dims=1)
    #Prod_Tk=reshape(Prod_Tk,N,N)'
    #temp=prod(temp,dims=1)
    #temp=reshape(temp,N,N)
    #Sum_Tk=sum(temp,dims=1)'
    #temp=sum(temp,dims=1)
    #boundary_func=vcat(map(x->range^2-x^2, x))
    #func_value=prod(temp)*prod(boundary_func)
    return func_value,Sum_Tk,Prod_Tk,boundary_func
end

##切比雪夫多项式二阶导
function d2chebyshev(x,k)
    chebyshev_Tk=Jacobi.chebyshev(x, k)
    chebyshev_Uk=Jacobi.chebyshev2(x, k)
    if x!=1&-1
        d2_Tk=k*((k+1)*chebyshev_Tk-chebyshev_Uk)/(x^2-1)
    end
    if x==1
        d2_Tk=(k^4-k^2)/3
    end
    if x==-1
        d2_Tk=(-1)^k*(k^4-k^2)/3
    end
    return d2_Tk
end

##基函数laplace
function laplace_psi(N,d,x,k,range)
    Sum_laplace=0
    Prod_laplace=0
    basisfunc=basis_function(N,d,x,k,range)[1]
    Sum_Tk=basis_function(N,d,x,k,range)[2]
    Prod_Tk=basis_function(N,d,x,k,range)[3]
        for dimension=1:d
            for N_x=1:N
                for i=1:N-1
                #################d2T####################
                d2T_k=d2chebyshev(x[N_x,dimension], k[i,dimension])
                T_k=Jacobi.chebyshev(x[N_x,dimension], k[i,dimension])
                Prod_laplace=basisfunc*d2T_k*Prod_Tk[i,N_x]/(Sum_Tk[i]*T_k)
                Sum_laplace+=Prod_laplace
                #################dTdT####################
                dT_ki=Jacobi.dchebyshev(x[N_x,dimension], k[i,dimension])
                T_k_i=Jacobi.chebyshev(x[N_x,dimension], k[i,dimension])
                    for j=i+1:N
                    dT_kj=Jacobi.dchebyshev(x[N_x,dimension], k[j,dimension])
                    T_k_j=Jacobi.chebyshev(x[N_x,dimension], k[j,dimension])
                    Prod_laplace=basisfunc*dT_kj*dT_ki*Prod_Tk[i,N_x]*Prod_Tk[j,N_x]/(Sum_Tk[i]*Sum_Tk[j]*T_k_i*T_k_j)
                    Sum_laplace+=2*Prod_laplace
                    end
                #################dTdx####################
                Prod_laplace=-2*x[N_x,dimension]*basisfunc*dT_ki*Prod_Tk[i,N_x]/(Sum_Tk[i]*(range^2-x[N_x,dimension]^2)*T_k_i)
                Sum_laplace+=2*Prod_laplace
                end
            T_k=Jacobi.chebyshev(x[N_x,dimension], k[N,dimension])
            d2T_k=d2chebyshev(x[N_x,dimension], k[N,dimension])
            Prod_laplace=basisfunc*d2T_k*Prod_Tk[N,N_x]/(Sum_Tk[N]*T_k)
            Sum_laplace+=Prod_laplace
            T_k_i=Jacobi.chebyshev(x[N_x,dimension], k[N,dimension])
            dT_ki=Jacobi.dchebyshev(x[N_x,dimension], k[N,dimension])
            Prod_laplace=-2*basisfunc*dT_ki*x[N_x,dimension]*Prod_Tk[N,N_x]/(Sum_Tk[N]*(range^2-x[N_x,dimension]^2)*T_k_i)
            Sum_laplace+=2*Prod_laplace
            #################d2x####################
            Sum_laplace+=-2*basisfunc/(range^2-x[N_x,dimension]^2)
            end
        end
    return Sum_laplace
end

function V_psi(N,d,x,k,range)
    Sum_V=0
    Sum_r=sum(x.^2)
        Sum_V=basis_function(N,d,x,k,range)[1]*Sum_r
    return Sum_V
end

function H_psi(N,d,x,k,range)
       Sum_H_psi=-1/2*laplace_psi(N,d,x,k,range)+1/2*V_psi(N,d,x,k,range)
    return Sum_H_psi
end

#以M=3,N=2为例生成 
#3  0
#2  1
function nsumk(N,M)
    k=collect(partitions(M+N,N))
    k=vcat(map(e->collect(e)', k)...)
    k=k-ones(Int64,size(k))
    return k
end

function find_k(d,N,M)
    k_all=zeros(Int64,1,N*d)
    for i=1:M
        dividers=nsumk(N*d,i)
        k_all=vcat(k_all, dividers)
    end
    k_all=sort!(k_all,dims =2)
    k_all=k_all'
    number_of_k=size(k_all)[2]
    i=0
    k=collect(permutations(k_all[1+i*N*d:N*d+i*N*d]))
    k=vcat(map(e->collect(e)', k)...)
    for i=1:number_of_k-1
        k_temp=collect(permutations(k_all[1+i*N*d:N*d+i*N*d]))
        k_temp=vcat(map(e->collect(e)', k_temp)...)
        k=vcat(k, k_temp)
        k=unique(k,dims=1)
    end
    k=k'
    number_of_k=size(k)[2]
    k=reshape(k,number_of_k*N*d)
    k=reshape(k,N,d,number_of_k)
    for i=1:number_of_k
        k[:,:,i]=sortslices(k[:,:,i],dims=1)
    end
    k=unique(k,dims=3)  
   return k
end

function random_choose(Number_of_x,d,N,k,range)
    Number_of_k=size(k)[3]
    random_choose_X=rand(N,d,Number_of_x).*rand([-range,range],N,d,Number_of_x)
    random_choose_psi_x=zeros(Number_of_k,Number_of_x)
    random_choose_Hpsi_x=zeros(Number_of_k,Number_of_x)
    for i=1:Number_of_k
        for j=1:Number_of_x
            random_choose_psi_x[i,j]=basis_function(N,d,random_choose_X[:,:,j],k[:,:,i],range)[1]
            random_choose_Hpsi_x[i,j]=H_psi(N,d,random_choose_X[:,:,j],k[:,:,i],range)
        end
    end
    return random_choose_psi_x,random_choose_Hpsi_x
end

function least_squares(least_squares_x,least_squares_y)
    least_squares_C_k=inv(least_squares_x'*least_squares_x)*least_squares_x'*least_squares_y
    return least_squares_C_k
end

function Eigenvalue(random_x,coef,alpha,d,N,k,range,Number_of_x;eps_x = 10^(-20),eps_d=0.001,maxIterations = 100)
    record=spzeros(maxIterations)
    basisfunction=random_x[1]
    H_basisfunction=random_x[2]
    for t in 1:maxIterations
            #若更换选点，取消以下注释
            #step_interval=10
            #if mod(t,step_interval)==0
                #random_x=random_choose(Number_of_x,d,N,k,range)
                #basisfunction=random_x[1]
                #H_basisfunction=random_x[2]
            #end
            reilaygh_quotient=reilaygh_quotient_alpha(coef,basisfunction,H_basisfunction,N,d,range,Number_of_x)[1]
            H_wave=reilaygh_quotient_alpha(coef,basisfunction,H_basisfunction,N,d,range,Number_of_x)[2]
            coef=reilaygh_quotient_alpha(coef,basisfunction,H_basisfunction,N,d,range,Number_of_x)[3]
            p=reilaygh_quotient_alpha(coef,basisfunction,H_basisfunction,N,d,range,Number_of_x)[4]
            coef_H_wave=reilaygh_quotient_alpha(coef,basisfunction,H_basisfunction,N,d,range,Number_of_x)[5]
            record[t]=reilaygh_quotient[1]
        if t>=3
            if abs(record[t]-record[t-1])+abs(record[t-1]-record[t-2])<=2*eps_d
                println("Convergence is reached after  ",t,"  iterations.")
                return coef,reilaygh_quotient,record
            end
        end
        alpha=steep(p,coef,reilaygh_quotient,coef_H_wave,basisfunction,N,d,range,H_basisfunction,Number_of_x,0,0.1,0.01,100)[1]
        println("============================步长为==========================",alpha)
        coef = coef-alpha*p
        println("At step ", t, " and reilaygh quotient = ", reilaygh_quotient)
    end
    println("Warning:",maxIterations,"  iterations have been exceeded")
    return coef,reilaygh_quotient,record
end

function reilaygh_quotient_alpha(coef_temp,basisfunction,H_basisfunction,N,d,range,Number_of_x)
        wave_temp=coef_temp*basisfunction
        sum_2_wave_temp=wave_temp*wave_temp'
        integal_wave_temp=(2*range)^(N*d)*sum_2_wave_temp/Number_of_x
        coef_temp=coef_temp/sqrt(integal_wave_temp[1])
        wave_temp=wave_temp/sqrt(integal_wave_temp[1])
        sum_2_wave_temp=sum_2_wave_temp/integal_wave_temp[1]
        H_wave_temp=coef_temp*H_basisfunction
        sum_wave_H_wave_temp=wave_temp*H_wave_temp'
        reilaygh_quotient_temp=sum_wave_H_wave_temp/sum_2_wave_temp
        coef_H_wave_temp=least_squares(basisfunction',H_wave_temp')
        coef_H_wave_temp=coef_H_wave_temp'
        p=coef_H_wave_temp-reilaygh_quotient_temp*coef_temp
    return reilaygh_quotient_temp,H_wave_temp,coef_temp,p,coef_H_wave_temp
end

function steep(p,coef,reilaygh_quotient,coef_H_wave,basisfunction,N,d,range,H_basisfunction,Number_of_x,a,b,eps_alpha,max)
        a_n=a
        b_n=b
        lambda=a+0.382*(b-a)
        mu=a+0.618*(b-a)
        phi_a=reilaygh_quotient_alpha(coef-a*p,basisfunction,H_basisfunction,N,d,range,Number_of_x)[1]
        phi_b=reilaygh_quotient_alpha(coef-b*p,basisfunction,H_basisfunction,N,d,range,Number_of_x)[1]
        phi_lambda=reilaygh_quotient_alpha(coef-lambda*p,basisfunction,H_basisfunction,N,d,range,Number_of_x)[1]
        phi_mu=reilaygh_quotient_alpha(coef-mu*p,basisfunction,H_basisfunction,N,d,range,Number_of_x)[1]
        phi=findmin([phi_a[1],phi_lambda[1],phi_mu[1],phi_b[1]])[1]
        t=findmin([phi_a[1],phi_lambda[1],phi_mu[1],phi_b[1]])[2]
        phi_max=findmax([phi_a[1],phi_lambda[1],phi_mu[1],phi_b[1]])[1]
    for i=1:max
        if t<3
            if b_n-a_n<eps_alpha
                return lambda
            end
            b_n=mu
            mu=lambda
            phi_b=phi_mu
            phi_mu=phi_lambda
            lambda=a_n+0.382*(b_n-a_n)
            phi_lambda=reilaygh_quotient_alpha(coef-lambda*p,basisfunction,H_basisfunction,N,d,range,Number_of_x)[1]
            phi=findmin([phi_a[1],phi_lambda[1],phi_mu[1],phi_b[1]])[1]
            t=findmin([phi_a[1],phi_lambda[1],phi_mu[1],phi_b[1]])[2]
            phi_max=findmax([phi_a[1],phi_lambda[1],phi_mu[1],phi_b[1]])[1]
        end
        if t>=3
           if b_n-a_n<eps_alpha
                return mu
           end
                a_n=lambda
                lambda=mu
                phi_a=phi_lambda
                phi_lambda=phi_mu
                mu=a_n+0.618*(b_n-a_n)
                phi_mu=reilaygh_quotient_alpha(coef-mu*p,basisfunction,H_basisfunction,N,d,range,Number_of_x)[1]
                phi=findmin([phi_a[1],phi_lambda[1],phi_mu[1],phi_b[1]])[1]
                t=findmin([phi_a[1],phi_lambda[1],phi_mu[1],phi_b[1]])[2]
                phi_max=findmax([phi_a[1],phi_lambda[1],phi_mu[1],phi_b[1]])[1]
        end
    end
    println("Warning:",max,"  iterations have been exceeded")
    return lambda
end

##########################
#求 d=1 N=3 M=2:20 reilaygh quotient变化
#运行时间较长
##########################
d=1
N=3
M=2
k=find_k(d,N,M)
range=5
Number_k=size(k)[3]
c_k=ones(1,Number_k)
eps_x = 10^(-6)
eps_d=0.001
alpha=0.001
Number_of_x=40000
println("=================================M=",M,"==============================================")
println("========================选点并计算基函数值时间如下","==============================================")
@time random_x=random_choose(Number_of_x,d,N,k,range)
@time df1=Eigenvalue(random_x,c_k,alpha,d,N,k,range,Number_of_x;eps_x= 10^(-30),eps_d=10^(-7),maxIterations = 100000)[3]    
println("========================梯度下降时间如上","==============================================")
for M in 4:2:20
    println("=================================M=",M,"==============================================")
    k=find_k(d,N,M)
    Number_k=size(k)[3]
    c_k=ones(1,Number_k)
    println("========================选点并计算基函数值时间如下","==============================================")
    @time random_x=random_choose(Number_of_x,d,N,k,range)
    @time df2=Eigenvalue(random_x,c_k,alpha,d,N,k,range,Number_of_x;eps_x= 10^(-30),eps_d=10^(-7),maxIterations = 100000)[3]  
    println("========================梯度下降时间如上","==============================================")
    record=hcat(df1, df2)
    df1=record
end

#画时间图
#没写函数，手动改时间。。
x=collect(2:2:20)
time=[0.293095+0.429859,1.313301+1.085863,2.614691+5.175826,4.714854+14.778792,7.623474+30.454066,12.123610+59.423014,18.228154+100.376995,25.576048+176.916978,36.163362+342.453232,43.929436+1113.821112 ]
Plots.plot(x,time,linewidth=2,title="d=1,N=3",xlabel="M",ylabel="time/second",label="",color="orange")
scatter!(x,time,label="",color="orange")

##改true_value 
true_value=1.5
record=df1
i=1
reilaygh_quotient_array1=findall(!iszero, record[:,i])
step=size(reilaygh_quotient_array1)
reilaygh_quotient_array1=Array(record[:,i][1:step[1]])-true_value*ones(step[1],1)
reilaygh_quotient_array1_ln=map((x) -> log(abs(x)), reilaygh_quotient_array1)

i=2
reilaygh_quotient_array2=findall(!iszero, record[:,i])
step=size(reilaygh_quotient_array2)
reilaygh_quotient_array2=Array(record[:,i][1:step[1]])-true_value*ones(step[1],1)
reilaygh_quotient_array2_ln=map((x) -> log(abs(x)), reilaygh_quotient_array2)

i=3
reilaygh_quotient_array3=findall(!iszero, record[:,i])
step=size(reilaygh_quotient_array3)
reilaygh_quotient_array3=Array(record[:,i][1:step[1]])-true_value*ones(step[1],1)
reilaygh_quotient_array3_ln=map((x) -> log(abs(x)), reilaygh_quotient_array3)

i=4
reilaygh_quotient_array4=findall(!iszero, record[:,i])
step=size(reilaygh_quotient_array4)
reilaygh_quotient_array4=Array(record[:,i][1:step[1]])-true_value*ones(step[1],1)
reilaygh_quotient_array4_ln=map((x) -> log(abs(x)), reilaygh_quotient_array4)

i=5
reilaygh_quotient_array5=findall(!iszero, record[:,i])
step=size(reilaygh_quotient_array5)
reilaygh_quotient_array5=Array(record[:,i][1:step[1]])-true_value*ones(step[1],1)
reilaygh_quotient_array5_ln=map((x) -> log(abs(x)), reilaygh_quotient_array5)

i=6
reilaygh_quotient_array6=findall(!iszero, record[:,i])
step=size(reilaygh_quotient_array6)
reilaygh_quotient_array6=Array(record[:,i][1:step[1]])-true_value*ones(step[1],1)
reilaygh_quotient_array6_ln=map((x) -> log(abs(x)), reilaygh_quotient_array6)

i=7
reilaygh_quotient_array7=findall(!iszero, record[:,i])
step=size(reilaygh_quotient_array7)
reilaygh_quotient_array7=Array(record[:,i][1:step[1]])-true_value*ones(step[1],1)
reilaygh_quotient_array7_ln=map((x) -> log(abs(x)), reilaygh_quotient_array7)

i=8
reilaygh_quotient_array8=findall(!iszero, record[:,i])
step=size(reilaygh_quotient_array8)
reilaygh_quotient_array8=Array(record[:,i][1:step[1]])-true_value*ones(step[1],1)
reilaygh_quotient_array8_ln=map((x) -> log(abs(x)), reilaygh_quotient_array8)

i=9
reilaygh_quotient_array9=findall(!iszero, record[:,i])
step=size(reilaygh_quotient_array9)
reilaygh_quotient_array9=Array(record[:,i][1:step[1]])-true_value*ones(step[1],1)
reilaygh_quotient_array9_ln=map((x) -> log(abs(x)), reilaygh_quotient_array9)

i=10
reilaygh_quotient_array10=findall(!iszero, record[:,i])
step=size(reilaygh_quotient_array10)
reilaygh_quotient_array10=Array(record[:,i][1:step[1]])-true_value*ones(step[1],1)
reilaygh_quotient_array10_ln=map((x) -> log(abs(x)), reilaygh_quotient_array10)


Plots.plot(reilaygh_quotient_array1,linewidth=2,title="d=1,N=3",xlabel="steps",ylabel="err",label="M=2",color="orange")
#scatter!(reilaygh_quotient_array1,label="",color="orange")

plot!(reilaygh_quotient_array2,linewidth=2,label="M=4",color="darkorange2")
#scatter!(reilaygh_quotient_array2,label="",color="darkorange2")

plot!(reilaygh_quotient_array3,linewidth=2,label="M=6",color="orange3")
#scatter!(reilaygh_quotient_array3,label="",color="orange3")

plot!(reilaygh_quotient_array4,linewidth=2,label="M=8",color="orangered")
#scatter!(reilaygh_quotient_array4,label="",color="orangered")

plot!(reilaygh_quotient_array5,linewidth=2,label="M=10",color="firebrick")
#scatter!(reilaygh_quotient_array5,label="",color="firebrick")

plot!(reilaygh_quotient_array6,linewidth=2,label="M=12",color="red4")
#scatter!(reilaygh_quotient_array6,label="",color="red4")

plot!(reilaygh_quotient_array7,linewidth=2,label="M=14",color="lightslateblue")
#scatter!(reilaygh_quotient_array7,label="",color="lightslateblue")

plot!(reilaygh_quotient_array8,linewidth=2,label="M=16",color="indigo")
#scatter!(reilaygh_quotient_array8,label="",color="indigo")

plot!(reilaygh_quotient_array9,linewidth=2,label="M=18",color="darkslategrey")
#scatter!(reilaygh_quotient_array9,label="",color="darkslategrey")

plot!(reilaygh_quotient_array10,linewidth=2,label="M=20",color="navy")
#scatter!(reilaygh_quotient_array10,label="",color="navy")

#plot(zeros(1:30),linewidth=1,ls=:dash,color="gray",label="")

Plots.plot(reilaygh_quotient_array1_ln,linewidth=2,title="d=1,N=3",xlabel="steps",ylabel="ln err",label="M=2",color="orange")
#scatter!(reilaygh_quotient_array1_ln,label="",color="orange")

plot!(reilaygh_quotient_array2_ln,linewidth=2,label="M=4",color="darkorange2")
#scatter!(reilaygh_quotient_array2_ln,label="",color="darkorange2")

plot!(reilaygh_quotient_array3_ln,linewidth=2,label="M=6",color="orange3")
#scatter!(reilaygh_quotient_array3_ln,label="",color="orange3")

plot!(reilaygh_quotient_array4_ln,linewidth=2,label="M=8",color="orangered")
#scatter!(reilaygh_quotient_array4_ln,label="",color="orangered")

plot!(reilaygh_quotient_array5_ln,linewidth=2,label="M=10",color="firebrick")
#scatter!(reilaygh_quotient_array5_ln,label="",color="firebrick")

plot!(reilaygh_quotient_array6_ln,linewidth=2,label="M=12",color="red4")
#scatter!(reilaygh_quotient_array6_ln,label="",color="red4")

plot!(reilaygh_quotient_array7_ln,linewidth=2,label="M=14",color="lightslateblue",title="d=1,N=2",xlabel="steps",ylabel="ln reilaygh quotient err")
#scatter!(reilaygh_quotient_array7_ln,label="",color="lightslateblue")

plot!(reilaygh_quotient_array8_ln,linewidth=2,label="M=16",color="indigo")
#scatter!(reilaygh_quotient_array8_ln,label="",color="indigo")

plot!(reilaygh_quotient_array9_ln,linewidth=2,label="M=18",color="darkslategrey")
#scatter!(reilaygh_quotient_array9_ln,label="",color="darkslategrey")

plot!(reilaygh_quotient_array10_ln,linewidth=2,label="M=20",color="navy")
#scatter!(reilaygh_quotient_array10_ln,label="",color="navy")

Plots.savefig("d=1N=3M=2:20.png")

err=[reilaygh_quotient_array1[size(reilaygh_quotient_array1)[1]][1],
    reilaygh_quotient_array2[size(reilaygh_quotient_array2)[1]],
    reilaygh_quotient_array3[size(reilaygh_quotient_array3)[1]],
    reilaygh_quotient_array4[size(reilaygh_quotient_array4)[1]],
    reilaygh_quotient_array5[size(reilaygh_quotient_array5)[1]],
    reilaygh_quotient_array6[size(reilaygh_quotient_array6)[1]],
    reilaygh_quotient_array7[size(reilaygh_quotient_array7)[1]],
    reilaygh_quotient_array8[size(reilaygh_quotient_array8)[1]],
    reilaygh_quotient_array9[size(reilaygh_quotient_array9)[1]],
    reilaygh_quotient_array10[size(reilaygh_quotient_array10)[1]]]
ln_err=map(x->log(abs(x)),err)

x=collect(2:2:20)
y=collect(2:-1:-7)
Plots.plot(x,ln_err,linewidth=2,label="",color="orange",xlabel="M",ylabel="ln err",title="d=1,N=2")
plot!(x,y,linewidth=1,ls=:dash,color="gray",label="")


function Eigenvalue1(random_x,coef,alpha,d,N,k,range,Number_of_x;eps_x = 10^(-20),eps_d=0.001,maxIterations = 100)
    record=spzeros(maxIterations)
    basisfunction=random_x[1]
    H_basisfunction=random_x[2]
    for t in 1:maxIterations
            #若更换选点，取消以下注释
            #step_interval=10
            #if mod(t,step_interval)==0
                #random_x=random_choose(Number_of_x,d,N,k,range)
                #basisfunction=random_x[1]
                #H_basisfunction=random_x[2]
            #end
            reilaygh_quotient=reilaygh_quotient_alpha(coef,basisfunction,H_basisfunction,N,d,range,Number_of_x)[1]
            H_wave=reilaygh_quotient_alpha(coef,basisfunction,H_basisfunction,N,d,range,Number_of_x)[2]
            coef=reilaygh_quotient_alpha(coef,basisfunction,H_basisfunction,N,d,range,Number_of_x)[3]
            p=reilaygh_quotient_alpha(coef,basisfunction,H_basisfunction,N,d,range,Number_of_x)[4]
            coef_H_wave=reilaygh_quotient_alpha(coef,basisfunction,H_basisfunction,N,d,range,Number_of_x)[5]
            record[t]=reilaygh_quotient[1]
        if t>=3
            if abs(record[t]-record[t-1])+abs(record[t-1]-record[t-2])<=2*eps_d
                println("Convergence is reached after  ",t,"  iterations.")
                return coef,reilaygh_quotient,record
            end
        end
        alpha=steep(p,coef,reilaygh_quotient,coef_H_wave,basisfunction,N,d,range,H_basisfunction,Number_of_x,0,0.1,0.01,100)[1]
        println("============================步长为==========================",alpha)
        coef = coef-alpha*p
        println("At step ", t, " and reilaygh quotient = ", reilaygh_quotient)
    end
    println("Warning:",maxIterations,"  iterations have been exceeded")
    return coef,reilaygh_quotient,record
end

d=1
N=4
M=20
@time k=find_k(d,N,M)
Number_of_x=40000
range=5
@time random_x=random_choose(Number_of_x,d,N,k,range)

Number_k=size(k)[3]
coef=ones(1,Number_k)
eps_x = 10^(-6)
eps_d=0.001
alpha=0.001
@time A=Eigenvalue1(random_x,coef,alpha,d,N,k,range,Number_of_x;eps_x= 10^(-30),eps_d=10^(-7),maxIterations = 100000)

basisfunction=random_x[1]
H_basisfunction=random_x[2]
H_wave=reilaygh_quotient_alpha(coef,basisfunction,H_basisfunction,N,d,range,Number_of_x)[2]
coef=reilaygh_quotient_alpha(coef,basisfunction,H_basisfunction,N,d,range,Number_of_x)[3]
p=reilaygh_quotient_alpha(coef,basisfunction,H_basisfunction,N,d,range,Number_of_x)[4]
coef_H_wave=reilaygh_quotient_alpha(coef,basisfunction,H_basisfunction,N,d,range,Number_of_x)[5]
reilaygh_quotient=reilaygh_quotient_alpha(coef,basisfunction,H_basisfunction,N,d,range,Number_of_x)[1]

coef=A[1]
record_1=spzeros(10000)
t=1
for alpha=-0.2:0.01:0.2
    coef=A[1]
    coef=coef-alpha*p
    f_x=reilaygh_quotient_alpha(coef,basisfunction,H_basisfunction,N,d,range,Number_of_x)[1]
    record_1[t]=f_x[1]
    t=t+1
end

coef=A[1]
coef=reilaygh_quotient_alpha(coef,basisfunction,H_basisfunction,N,d,range,Number_of_x)[3]
p=reilaygh_quotient_alpha(coef,basisfunction,H_basisfunction,N,d,range,Number_of_x)[4]
record_1=spzeros(10000)
coef=coef-0.01*p
coef=reilaygh_quotient_alpha(coef,basisfunction,H_basisfunction,N,d,range,Number_of_x)[3]
p=reilaygh_quotient_alpha(coef,basisfunction,H_basisfunction,N,d,range,Number_of_x)[4]
t=1
for alpha=-0.2:0.01:0.2
    coef=reilaygh_quotient_alpha(coef,basisfunction,H_basisfunction,N,d,range,Number_of_x)[3]
    coef=coef-alpha*p
    f_x=reilaygh_quotient_alpha(coef,basisfunction,H_basisfunction,N,d,range,Number_of_x)[1]
    record_1[t]=f_x[1]
    t=t+1
end

reilaygh_quotient_array=findall(!iszero, record_1)
step=size(reilaygh_quotient_array)
reilaygh_quotient_array=Array(record_1[1:step[1]])
findmin(reilaygh_quotient_array)

x=collect(-0.2:0.01:0.2)
Plots.plot(x,reilaygh_quotient_array,linewidth=2,ylabel="f(x)",label="",xlabel="alpha")
