import math

def generate_params_mix_net():
    print("##### Decryption Proof for Mix-Net Based Voting #####")
    N = 4096
    tau = 4096
    q = 4611686018326724609
    p = 2 
    lambda_ = 128
    k = 3
    n = 1
    v = k + 1  # in AS=T equation, s_i \in R_q^v
    ntil = 130
    B = 1
    eta = 36

    sigma_A = 0.675 * (B*k*N + (2*N+1)*N)*tau
    Bound_A = math.sqrt(2*v*N) * sigma_A
    print("sigma_A: ", sigma_A, math.log2(sigma_A))
    print("Bound_A bits: ", Bound_A, math.log2(Bound_A))

    sigma_C = 11 * eta * B * math.sqrt(k*N)
    Bound_C = 2 * sigma_C * math.sqrt(N)
    print("sigma_C: ", sigma_C, math.log2(sigma_C))
    print("Bound_C: ", Bound_C, math.log2(Bound_C))

    KB = 8 * 1024
    q = 2**64 # using 64 bits to store coefficients

    message_size = N*math.log2(p) / KB
    print("m_i: ", message_size, "KB")

    ciphertext_size = 2*N*math.log2(q) / KB
    print("(u_i,v_i): ", ciphertext_size, "KB")

    commitment_size = 2*n*N*math.log2(q) / KB
    print("d_i: ", commitment_size, "KB")

    lin_proof_size = (2*(k- n)*N*math.log2(6*sigma_C) / KB) 
    print("pi_lin: ", lin_proof_size, "KB")

    shortness_proof_size = ((k+1)*ntil*N*math.log2(6*sigma_A) / KB) / tau
    print("pi_a: ", shortness_proof_size, "KB") 

    dec_proof_size = shortness_proof_size  + (lin_proof_size + commitment_size) 
    print("pi_dec: ", dec_proof_size , "KB") 



def generate_params_homomorphic():
    print("##### Decryption Proof for Homomorphic Tallying Voting #####")
    N = 4096
    tau = 1000000
    q = 4611686018326724609
    p = 10000009
    lambda_ = 128
    k = 3
    n = 1
    v = k + 1  # in AS=T equation, s_i \in R_q^v
    ntil = 130
    B = 1
    eta = 36

    Bound_A = q / (4*p)
    sigma_A = Bound_A / math.sqrt(2*v*N)
    print("sigma_A: ", sigma_A, math.log2(sigma_A))
    print("Bound_A bits: ", Bound_A, math.log2(Bound_A))

    sigma_C = 11 * eta * B * math.sqrt(k*N)
    Bound_C = 2 * sigma_C * math.sqrt(N)
    print("sigma_C: ", sigma_C, math.log2(sigma_C))
    print("Bound_C: ", Bound_C, math.log2(Bound_C))

    KB = 8 * 1024
    q = 2**64 # using 64 bits to store coefficients

    message_size = N*math.log2(p) / KB
    print("m_i: ", message_size, "KB")

    ciphertext_size = 2*N*math.log2(q) / KB
    print("(u_i,v_i): ", ciphertext_size, "KB")

    commitment_size = 2*n*N*math.log2(q) / KB
    print("d_i: ", commitment_size, "KB")

    lin_proof_size = (2*(k- n)*N*math.log2(6*sigma_C) / KB) 
    print("pi_lin: ", lin_proof_size, "KB")

    shortness_proof_size = ((k+1)*ntil*N*math.log2(6*sigma_A) / KB)
    print("pi_a: ", shortness_proof_size, "KB") 

    dec_proof_size = shortness_proof_size  + (lin_proof_size + commitment_size) 
    print("pi_dec: ", dec_proof_size / tau, "KB") 



generate_params_mix_net()
generate_params_homomorphic()