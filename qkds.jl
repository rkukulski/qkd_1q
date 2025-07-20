
function SixStateProtocol()
    #SSP/Sarg
    #https://doi.org/10.1103/PhysRevLett.81.3018
    name = "SixStateProtocol"
    qa = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6]
    list_rho = [[1, 0], [0, 1], [1/sqrt(2), 1/sqrt(2)], [1/sqrt(2), -1/sqrt(2)], [1/sqrt(2), im/sqrt(2)], [1/sqrt(2), -im/sqrt(2)]]
    povm = [[0, 1]/sqrt(3), [1, 0]/sqrt(3), [1, -1]/sqrt(6), [1, 1]/sqrt(6), [1, -im]/sqrt(6), [1, im]/sqrt(6)]
    f = [[1/3, 1/3, 1/3], [[1,2], [3,4], [5,6]]]
    return [qa, [x*x' for x in list_rho], [x*x' for x in povm], f], name
end
function MSZProtocol()
    #Nie są stany podane, tylko, że nieortogonalne
    #Zaimplememtowałem to na podstawie fig. 2
    #Raczej nie zadziala
    #https://doi.org/10.1016/0030-4018(95)00688-5
    name = "MSZProtocol"
    qa = [1/4, 1/4, 1/4, 1/4]
    list_rho = [[1, 0], [-1, 0], [0, 1], [0, -1]]
    povm = [[-1, 0]/sqrt(2), [1, 0]/sqrt(2), [0, -1]/sqrt(2), [0, 1]/sqrt(2)]
    f = [[1/2, 1/2],[[1,2], [3,4]]]
    return [qa, [x*x' for x in list_rho], [x*x' for x in povm], f], name
end
#S13 ciekawe ale się nie da zaimplementować
function KMB09_N4()
    #To nie jest do konca ich propozycja, ale niech bedzie
    name = "KMB09_N4"
    qa = [1/4, 1/4, 1/4, 1/4]
    list_rho = [
        [1/2, 1/2, 1/2, 1/2],
        [1/2, 1/2, -1/2, -1/2],
        [1/2, -1/2, 1/2, -1/2],
        [1/2, -1/2, -1/2, 1/2]
    ]
  
    povm = [
        [1/2, 1/2, -1/2, -1/2],    
        [1/2, -1/2, 1/2, -1/2],   
        [1/2, -1/2, -1/2, 1/2],   
        [1/2, 1/2, 1/2, 1/2]       
    ]
    f = [[1.0],[[1,2,3,4]]] 
    return [qa, [x*x' for x in list_rho], [x*x' for x in povm], f], name
end

#W pewnym artykule taki asymetryczny sposob pomiaru jest tez stosowany niby w Twin Field
function T12(px::Float64 = 1/16)
    name = "T12"
#https://doi.org/10.1364/OE.21.024550
#Wybor bazy nie jest tak samo prawdopodobny jak w BB84
    pz = 1 - px
    qa = [
        pz * 0.5, 
        pz * 0.5, 
        px * 0.5, 
        px * 0.5  
    ]

    list_rho = [
        [1, 0],             
        [0, 1],              
        [1/sqrt(2), 1/sqrt(2)], 
        [1/sqrt(2), -1/sqrt(2)]  
    ]

  
    povm = [
        [0, 1]/sqrt(2),    
        [1, 0]/sqrt(2),     
        [1/2, -1/2],        
        [1/2, 1/2]           
    ]
    f = [[pz, px], [[1,2], [3,4]]]

    return [qa, [x*x' for x in list_rho], [x*x' for x in povm], f], name
end
function Chau15()
    #Inny kompletnie sposob kodowania dla n wymiarów, tutaj pewne uproszczenie dla 4 wymiarów
    #https://doi.org/10.1103/PhysRevA.92.062324
    name = "Chau15"
    qa = [1/4, 1/4, 1/4, 1/4]

    list_rho = [
        [1/sqrt(2), 1/sqrt(2), 0, 0],         # (|0⟩ + |1⟩)/√2
        [1/sqrt(2), -1/sqrt(2), 0, 0],        # (|0⟩ - |1⟩)/√2
        [0, 0, 1/sqrt(2), 1/sqrt(2)],         # (|2⟩ + |3⟩)/√2
        [0, 0, 1/sqrt(2), -1/sqrt(2)]         # (|2⟩ - |3⟩)/√2
    ]

    povm = [
        [1/sqrt(2), -1/sqrt(2), 0, 0],        
        [0, 0, 1/sqrt(2), 1/sqrt(2)],         
        [0, 0, 1/sqrt(2), -1/sqrt(2)],        
        [1/sqrt(2), 1/sqrt(2), 0, 0]          
    ]

    f = [[1.0],[[1,2,3,4]]]

    return [qa, [x*x' for x in list_rho], [x*x' for x in povm], f], name
end


#Dla k=3 jest R04 
#k = 3,5,7 a in [0,2pi)]
#https://journals.aps.org/pra/abstract/10.1103/PhysRevA.89.042330 to wydaje sie podobne?
function K_State_P1(K_val::Int, a_val::Float64)
    name = "K_State_P1$(K_val)_a$(round(a_val, digits=3))"
    #Nie wiem jak mozna zrobic P2, tam nie mamy juz povms a unambigious measurement
    #1/K prob.
    qa = fill(1.0/K_val, K_val)

    list_rho = Vector{Vector{ComplexF64}}(undef, K_val)
    for k = 0:K_val-1
        theta_k = (2.0 * pi * k) / K_val
        angle = a_val + theta_k
        # |ψk⟩ = cos(a+θk)|0⟩ + sin(a+θk)|1⟩
        list_rho[k+1] = [cos(angle), sin(angle)]
    end

    # POVMs  M_k = (2/K) |ψk_perp⟩⟨ψk_perp|.
    povm = Vector{Vector{ComplexF64}}(undef, K_val)
    for k = 0:K_val-1
        theta_k = (2.0 * pi * k) / K_val
        angle = a_val + theta_k
        # |ψk_perp⟩ = sin(a+θk)|0⟩ - cos(a+θk)|1⟩ (a state orthogonal to |ψk⟩)
        povm_vector = [sin(angle), -cos(angle)]
        povm[k+1] = povm_vector * sqrt(2.0/K_val)
    end

    f = [[1.0], [collect(1:K_val)]]

    return [qa, [x*x' for x in list_rho], [x*x' for x in povm], f], name
end

function ThreeStateQKD()
    #http://dx.doi.org/10.1080/09500340008244056
    #dla stanu 1. też poodaje p=2/3
    name = "Three-state QKD "
    qa = [1/3, 1/3, 1/3]

    list_rho = [
        [1.0 + 0.0im, 0.0 + 0.0im],                          
        [-0.5 + 0.0im, sqrt(3)/2 + 0.0im],                  
        [-0.5 + 0.0im, -sqrt(3)/2 + 0.0im]                  
    ]

    povm = [
        sqrt(2/3) * [sqrt(3)/2 + 0.0im, -0.5 + 0.0im],
        sqrt(2/3) * [0.0 + 0.0im, 1.0 + 0.0im],
        sqrt(2/3) * [-sqrt(3)/2 + 0.0im, -0.5 + 0.0im],

    ]

    f = [[1.0],   
         [[0,1,2]]] 

    return [qa, [x*x' for x in list_rho], [x*x' for x in povm], f], name
end

#https://link.springer.com/article/10.1007/s11128-020-02927-8
#chyba nie da sie tego zaimplementować
# function QutritSQKD()
#     name = "QutritSQKD"

#     omega = exp(im * 2 * pi / 3)

#     q0 = [1.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im]
#     q1 = [0.0 + 0.0im, 1.0 + 0.0im, 0.0 + 0.0im]
#     q2 = [0.0 + 0.0im, 0.0 + 0.0im, 1.0 + 0.0im]
#     q01  = (1/sqrt(3)) * [omega, 1.0 + 0.0im, 1.0 + 0.0im]
#     q11 = (1/sqrt(3)) * [1.0 + 0.0im, omega, 1.0 + 0.0im]
#     q21  = (1/sqrt(3)) * [1.0 + 0.0im, 1.0 + 0.0im, omega]
#     q02 = (1/sqrt(3)) * [1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im]
#     q12 = (1/sqrt(3)) * [1.0 + 0.0im, omega, -omega]
#     q22 = (1/sqrt(3)) * [1.0 + 0.0im, -omega, omega]
#     list_rho = [q0, q1, q2, q01, q11, q21,q02, q12, q22]

#     qa = [1.0 / 9 for _ in 1:length(list_rho)]

#     povm_states = [
#         [1.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im],
#         [0.0 + 0.0im, 1.0 + 0.0im, 0.0 + 0.0im],
#         [0.0 + 0.0im, 0.0 + 0.0im, 1.0 + 0.0im]  
#     ]

#     f_distribution = [1/2, 1/2]

#     f_points = [[1, 2, 3,4,5,6], [1,2,3,7,8,9]]
#     f = [f_distribution, f_points]

#     return [qa, [x*x' for x in list_rho], [x*x' for x in povm_states], f], name
# end