This folder includes two example implementations:
(1) Optimized implementation of SPNbox-16:
      # OPTI_SPNBOX16_SBOXGEN.C: optimized sbox generation of SPNbox-16
      # OPTI_SPNBOX16_BBI.C: optimized black-box implementation of SPNbox-16
      # OPTI_SPNBOX16_WBI.C: optimized white-box implementation of SPNbox-16
      # Makefile, Makefile.am, Makefile.in: compile related files

      # sbox.h: original sbox of SPNbox-16
      # invsbox.h: original inverse sbox of SPNbox-16
      # invsbox01.h ~ invsbox10.h:  perturbed dysfunctional inverse sbox of SPNbox-16 (for 10 users)
     
      # pertu_ciphertext01.txt ~ pertu_ciphertext10.txt: randomly generated perturbation ciphertexts (for 10 users)
      # pertu_plaintext01.txt ~ pertu_plaintext10.txt: plaintexts of the perturbation ciphertexts (for 10 users)
      # newpertu_plaintext01.txt ~ newpertu_plaintext10.txt: plaintexts of the perturbation ciphertexts under the perturbed white-box decryption programs (for 10 users) 
      # random_ciphertext.txt: a randomly generated ciphertext
      # random_plaintext.txt: plaintext of the random ciphertext
      # newrandom_plaintext: plaintext of the random ciphertext under a perturbed white-box decryption program
      # testtraitor_plaintext01.txt ~ testtraitor_plaintext10.txt: test output of perturbation ciphertexts under the pirate decoder (we assume the 5th user is the traitor)

(2) Optimized implementation of WEM-16:
      # OPTI_WEM16_SBOXGEN.C: optimized sbox generation of WEM-16
      # OPTI_WEM16_WBI.C: optimized white-box implementation of WEM-16
      # Makefile, Makefile.am, Makefile.in: compile related files

      # sbox.h: original sbox of WEM-16
      # invsbox.h: original inverse sbox of WEM-16
      # newinvsbox01.h ~ newinvsbox10.h:  perturbed dysfunctional inverse sbox of WEM-16 (for 10 users)
     
      # pertu_ciphertext01.txt ~ pertu_ciphertext10.txt: randomly generated perturbation ciphertexts (for 10 users)
      # pertu_plaintext01.txt ~ pertu_plaintext10.txt: plaintexts of the perturbation ciphertexts (for 10 users)
      # newpertu_plaintext01.txt ~ newpertu_plaintext10.txt: plaintexts of the perturbation ciphertexts under the perturbed white-box decryption programs (for 10 users)
      # random_ciphertext.txt: a randomly generated ciphertext
      # random_plaintext.txt: plaintext of the random ciphertext
      # newrandom_plaintext: plaintext of the random ciphertext under a perturbed white-box decryption program
      # random_plaintext01.txt ~ random_plaintext10.txt : plaintext of the random ciphertext under perturbed white-box decryption programs (for 10 users)
      # testtraitor_plaintext01.txt ~ testtraitor_plaintext10.txt: test output of perturbation ciphertexts under the pirate decoder (we assume the 5th user is the traitor)


NOTE: These implementations need GIVARO (https://github.com/linbox-team/givaro), AES-NI, AVX2 and SSE2
ANY QUESTION PLEASE CONTACT: junl1212@163.com
