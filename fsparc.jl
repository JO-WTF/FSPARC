using JSON
using JuMP,Mosek,Ipopt,CPLEX


#################### Read data from data.json ####################
f=open("data.json","r")
data=readline(f)
data=JSON.parse(data)

P=data["P"]
I=data["I"]
J=data["J"]
bw=data["Bandwidth"]
ns=data["Noise"]
el=data["ell"]


probstatus=""

LB=0
UB=0
#################### Solve NLP with BONMIN to get UB ####################
nlp = Model(solver=IpoptSolver(print_level=0,max_iter=500000000,bound_relax_factor=1e-18))

@variable(nlp,nlp_p[1:I]>=0)
@variable(nlp,nlp_r[1:I]>=0)
@constraint(nlp,sum(nlp_p)<=P)
@NLconstraint(nlp,[i=1:I],nlp_r[i]<=bw[i]*log2(1+nlp_p[i]/ns[i]))

@objective(nlp,Max,sum(nlp_r))
tic()
nlpstatus=solve(nlp;suppress_warnings=true)
nlptime=toq();
milptime=0
if nlpstatus != :Infeasible

	#################### Get UB of the problem ####################
	pre_nlp_p_opt=getvalue(nlp_p)
	pre_nlp_r_opt=getvalue(nlp_r)
	pre_nlp_obj_opt=getobjectivevalue(nlp)
	if el=="NULL"
		sumz=0
		z=zeros(J)
		for j=1:J
			z[j]=2.71828^randn()
			sumz=sumz+z[j]
		end
		el=zeros(J)
		for j=1:J
			el[j]=z[j]/sumz*pre_nlp_obj_opt*DR
			# print(el[j],", ")
		end
		close(f)
	end

	print("The UB of the relaxed NLP is ", pre_nlp_obj_opt," and the DR is ", round(sum(el)/pre_nlp_obj_opt*100,2) ,"%\n")
	UB=pre_nlp_obj_opt
	if pre_nlp_obj_opt>=sum(el)
		
		#################### Pre-Processing 0-1 Mixed LP ####################
		mlp=Model(solver=CplexSolver(CPX_PARAM_MIPDISPLAY=0,CPX_PARAM_FPHEUR=1,CPX_PARAM_TILIM=5))
		@variable(mlp,mlp_x[1:I,1:J],Bin)
		@constraint(mlp,[i=1:I],sum(mlp_x[i,j] for j=1:J)<=1)
		@constraint(mlp,[j=1:J],sum(pre_nlp_r_opt[i]*mlp_x[i,j] for i=1:I)>=el[j])
		@objective(mlp,Max,sum(pre_nlp_r_opt[i]*mlp_x[i,j] for i=1:I, j=1:J))
		tic()
		mlpstatus=solve(mlp;suppress_warnings=false)
		mlptime=toq();
		mlpstatus=:Infeasible
		if mlpstatus==:Infeasible || mlpstatus==:UserLimit
			print("Unable to solve the problem with Pre-Processing method\n")
			#################### OA algorithm ####################
			tic()
			timelimit=120
			@printf("Time limit is: %.2fs\n",timelimit)
			milp=Model(solver=CplexSolver(CPX_PARAM_MIPDISPLAY=0,CPX_PARAM_FPHEUR=1,CPX_PARAM_TILIM=timelimit))
			@variable(milp,milp_x[1:I,1:J],Bin)
			@variable(milp,milp_p[1:I,1:J]>=0)
			@variable(milp,milp_r[1:I,1:J]>=0)

			@constraint(milp,sum(milp_p)<=P)
			@constraint(milp,[i=1:I],sum(milp_x[i,j] for j=1:J)<=1)
			@constraint(milp,[i=1:I,j=1:J],milp_p[i,j]<=P*milp_x[i,j])
			@constraint(milp,[j=1:J],sum(milp_r[i,j] for i=1:I)>=el[j])

			# warm-starting
			p_bar=fill(P,I,J)
			@constraint(milp,[i=1:I,j=1:J],milp_r[i,j]<= bw[i]/(log(2)*(ns[i]+p_bar[i,j]))*milp_p[i,j] + (bw[i]*log2(1+p_bar[i,j]/ns[i]) - bw[i]/(log(2)*(ns[i]+p_bar[i,j]))*p_bar[i,j])*milp_x[i,j])

			for i=1:I
				p_bar[i,:]=pre_nlp_p_opt[i]
			end

			@constraint(milp,[i=1:I,j=1:J],milp_r[i,j]<= bw[i]/(log(2)*(ns[i]+p_bar[i,j]))*milp_p[i,j] + (bw[i]*log2(1+p_bar[i,j]/ns[i]) - bw[i]/(log(2)*(ns[i]+p_bar[i,j]))*p_bar[i,j])*milp_x[i,j])

			for i=1:I
				p_bar[i,:]=sqrt(P*pre_nlp_p_opt[i])
			end
			@constraint(milp,[i=1:I,j=1:J],milp_r[i,j]<= bw[i]/(log(2)*(ns[i]+p_bar[i,j]))*milp_p[i,j] + (bw[i]*log2(1+p_bar[i,j]/ns[i]) - bw[i]/(log(2)*(ns[i]+p_bar[i,j]))*p_bar[i,j])*milp_x[i,j])


			@objective(milp,Max,sum(milp_r))
			print("Initial Formulation: ")
			# print(milp)
			milpstatus=solve(milp;suppress_warnings=false)
			timelimit=timelimit-toq()
			violated=true
			iteration=1;
			if milpstatus!=:Infeasible && milpstatus!=:InfeasibleOrUnbounded && milpstatus!=:UserLimit

				while violated==true && (sum(el[j] for j=1:J)-sum(getvalue(milp_r[i,j]) for i=1:I,j=1:J))/sum(el[j] for j=1:J)<=0.001

					#################### check violated constraints and add cutting planes ####################
					flag=0
					# print("x[i,j]","\t\t\t","p[i,j]","\t\t\t","r[i,j]\n")
					for i=1:I
						# for j=1:J
						# 	@printf("%d",getvalue(milp_x[i,j]))
						# 	print("\t")
						# end
						# print("\t")
						# for j=1:J
						# 	@printf("%.2f",getvalue(milp_p[i,j]))
						# 	print("\t")
						# end
						# print("\t")
						# for j=1:J
						# 	@printf("%.2f",getvalue(milp_r[i,j]))
						# 	print("\t")
						# end
						for j=1:J
							if getvalue(milp_x[i,j])>=0.5
								if getvalue(milp_r[i,j])-bw[i]*log2(1+getvalue(milp_p[i,j])/ns[i])>0.001
									p_bar[i,:]=ns[i]*(2^(getvalue(milp_r[i,j])/bw[i])-1)
									for j_prime=1:J
										@constraint(milp,milp_r[i,j_prime]<= bw[i]/(log(2)*(ns[i]+p_bar[i,j]))*milp_p[i,j_prime] + (bw[i]*log2(1+p_bar[i,j]/ns[i]) - bw[i]/(log(2)*(ns[i]+p_bar[i,j]))*p_bar[i,j])*milp_x[i,j_prime])
									end
									flag=1
									# @printf("r[%d,%d] is violated (%.5f > %.5f).",i,j,getvalue(milp_r[i,j]),bw[i]*log2(1+getvalue(milp_p[i,j])/ns[i]))
								end
							end
						end
						# println()
					end
					# print("current number of constraints:",MathProgBase.numconstr(milp),"\n")

					@printf("Current UB is %.6f.\n",sum(getvalue(milp_r)))
					UB=sum(getvalue(milp_r))
					if flag==1
						violated=true
						print("Iteration ",iteration,": ")
						iteration=iteration+1
						@printf("Time limit is: %.2fs\n",timelimit)
						if timelimit<0
							probstatus="Time-Limit"
							break
						end
						tic()
						setsolver(milp,CplexSolver(CPX_PARAM_MIPDISPLAY=0,CPX_PARAM_FPHEUR=1,CPX_PARAM_TILIM=timelimit))
						milpstatus=solve(milp;suppress_warnings=false)
						if milpstatus==:UserLimit||milpstatus==:Infeasible||milpstatus==:InfeasibleOrUnbounded
							timelimit=timelimit-toq()
							break
							probstatus="Time-Limit"
						end
						
						timelimit=timelimit-toq()

						################## Reallocate Power With X_ij Solved From MILP ####################

						renlp = Model(solver=IpoptSolver(print_level=0,max_iter=500000000,bound_relax_factor=1e-18))
						@variable(renlp,renlp_p[1:I]>=0)
						@variable(renlp,renlp_r[1:I]>=0)
						@constraint(renlp,sum(renlp_p)<=P)
						@NLconstraint(renlp,[i=1:I],renlp_r[i]<=bw[i]*log2(1+renlp_p[i]/ns[i]))
						channelsForUser=zeros(I,J)
						@objective(renlp,Max,sum(renlp_r))
						tic()
						renlpstatus=solve(renlp;suppress_warnings=true)
						toq()
						
						# for j=1:J
						# 	@printf("%.2f",sum(getvalue(renlp_r[i])*getvalue(milp_x[i,j]) for i=1:I))
						# 	print(" ")#
						# end
						# println()

						# flag=0
						# for j=1:J
						# 	if sum(getvalue(renlp_r[i])*getvalue(milp_x[i,j]) for i=1:I)<el[j]
						# 		flag=1
						# 	end
						# end
						# if flag==0
						# 	if sum(getvalue(renlp_r))>LB
						# 		LB=sum(getvalue(renlp_r))
						# 	end
						# 	if (UB-LB)/LB<0.001
						# 		violated=false
						# 		@printf("The problem is solved by renlp with total data rate: %.2f\n",sum(getvalue(renlp_r)))
						# 		probstatus="ReAlloc-Solved"
						# 	end
						# end


					else
						violated=false
						print("All demand satisfied, no cuts violated.\n")
						probstatus="OA-Solved"
					end

				end
				#end while violated==true
			else
				if milpstatus==:UserLimit
					print("Time limit reached.\n")
					probstatus="Time-Limit"
				else
					violated=false
					print("The problem is infeasible\n")
					probstatus="OA-Infeasible"
				end
			end
			#end if milpstatus!=:Infeasible && milpstatus!=:InfeasibleOrUnbounded
			

			milptime=120-timelimit

			@printf("MILP used %.2f seconds.\n",milptime)
		else
			if mlpstatus==:Optimal
				print("The problem has been solved by 0-1 Mixed LP.\n")
				probstatus="Pre-Processing"
			else
				println("The problem was not solved: ",mlpstatus,".\n")
			end
		end
		#end if mlpstatus==:Infeasible

	else
		print("The NLP problem is infeasible: Demand is too high.\n")
		probstatus="Deamnd-Too-High"
	end
	#end if pre_nlp_obj_opt>=sum(el)
else
	print("The NLP problem is infeasible. Something might have gone wrong.\n")
end
#end if nlpstatus != :Infeasible

# print(milp)

totaltime=round(nlptime+mlptime+milptime,2);

f=open("runningtime-single.txt","w")

write(f,"$totaltime\t$probstatus\n")

close(f)



