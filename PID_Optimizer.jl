using JuMP, Ipopt
function pid_optimizer(A,B,C,D,SollWert,max_overshoot,tf,Number_Of_Disc_Points,u_lower_bound,u_upper_bound)


    mod =Model(solver=IpoptSolver(max_iter=1000))
    System_dim = size(B)[1]-1
    A=reshape(A,System_dim+1,System_dim+1)'
    #Sytem variables
    @variable(mod,x[j=0:System_dim,i=0:Number_Of_Disc_Points])
    x0=zeros(System_dim+1,1)
    for ii =0:System_dim
        @constraint(mod,x[ii,0]==x0[ii+1])
    end
    #Integrator states
    @variable(mod,Integrator[i=0:Number_Of_Disc_Points])
    @constraint(mod,Integrator[0]==0)

    #Output
    @variable(mod,y[i=0:Number_Of_Disc_Points])
    @constraint(mod,y[0]==0)
    for ii =0:Number_Of_Disc_Points
        @constraint(mod,y[ii]<=max_overshoot)
    end
    #error
    @variable(mod,e[i=-1:Number_Of_Disc_Points])
    @constraint(mod,e[-1]==0)
    @constraint(mod,e[0]==SollWert-y[0])

    #PID Variables
    @variable(mod,K_P,start=50)
    @variable(mod,K_I,start=3)
    @variable(mod,K_D,start=3)
    #Formula for PID u[k] = K_P *e[k] + K_I sum_i=0^k(e[i])+ K_D * (e[k]-e[k-1])
    #SISO system
    @variable(mod,u_lower_bound<=u[i=0:Number_Of_Disc_Points]<=u_upper_bound)
    for ii =0:Number_Of_Disc_Points
        @constraint(mod,u[ii]==K_P*e[ii]+K_I*Integrator[ii]+K_D*(e[ii]-e[ii-1]))
    end

    #Kollokations bedingungen euler vorwärts
    delta_T=tf/Number_Of_Disc_Points
    for ii=0:Number_Of_Disc_Points-1
        for jj =0:System_dim
            temp=0
            for kk=0:System_dim
                temp = A[jj+1,kk+1]*x[kk,ii]+temp
            end
            @constraint(mod,x[jj,ii+1]-x[jj,ii] ==delta_T*( temp + B[jj+1]*u[ii]))
        end
        @constraint(mod,Integrator[ii+1]-Integrator[ii]==delta_T*e[ii])
    end
    for ii=0:Number_Of_Disc_Points
        temp2=0
        for jj=0:System_dim
            temp2 = C[jj+1]*x[jj,ii]+temp2
        end
        @constraint(mod,y[ii]==temp2+D[1]*u[ii])
        @constraint(mod,e[ii] == SollWert- y[ii])
    end



    @variable(mod,Cost)
    temp=0
    for ii=0:Number_Of_Disc_Points
        temp=temp+e[ii]*e[ii]
    end

    @constraint(mod,Cost== temp)
    @NLobjective(mod,:Min,Cost)
    status=solve(mod)
    getvalue(u)
    P=getvalue(K_P)
    I=getvalue(K_I)
    D=getvalue(K_D)
    e=getvalue(e)
    u=getvalue(u)
    y=getvalue(y)
    ret_y = zeros(Number_Of_Disc_Points+1)
    for ii=0:Number_Of_Disc_Points
        ret_y[ii+1] = y[ii]
    end
    x=getvalue(x)
    Integrator=getvalue(Integrator)
    println(P)
    println(I)
    println(D)
    #println(e[Number_Of_Disc_Points])
    t = linspace(0,tf,Number_Of_Disc_Points+1)
    #plot(t,y)#,color="red", linewidth=2.0, linestyle="--")
    #println(e)
    return status , [P,I,D],ret_y,t
end
