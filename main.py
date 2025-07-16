import random #simulating random selection process
from typing import Callable #type checking
import json #for caching
import os #for caching

class PRF:
  """
  A class that takes as input a security parameter and a key, generates and stores an n-bit safe prime p, where p = 2q + 1 and q is also prime,
  along with a generator g for the multiplicative group Z_p*.
  Here is an outline of the pseudorandom function:
    - It constructs a one way permutation using modular exponentiation of the generator g.
    - The hardcore predicate of the one way permutation which is the 
        xor of random (determined by the key) subset of bits of the input to the one way permutation.
    - A base pseudorandom generator (PRG) with expansion factor n + 1 is constructed using the one way permutation and its hardcore predicate.
    - Another PRG is constructed with expansion factor 2n using the base PRG.
    - Now, the pseudorandom function is constructed by using the PRG with expansion factor 2n .
        using the method described in Chapter 8 Section 8.5.
  """

  _parameter_cache = {}
  _cache_file = "dlp_parameters.json" # cache values of g, p, q for a given value of n for reuse

  @classmethod
  def _load_cache(cls):
    if os.path.exists(cls._cache_file):
      try:
        with open(cls._cache_file, 'r') as f:
          string_cache = json.load(f)
          return {int(k): tuple(v) for k, v in string_cache.items()}
      except:
        return {}
    return {}

  @classmethod
  def _save_cache(cls):
      with open(cls._cache_file, 'w') as f:
        # Need to convert values to lists since tuples aren't JSON serializable
        serializable_cache = {str(k): list(v) for k, v in cls._parameter_cache.items()}
        json.dump(serializable_cache, f)

  @staticmethod
  def _efficient_random_non_zero_bitstring(n: int):
    """
    Generates a random n-bit string that is not all zeroes.

    Args:
      n (int): Length of the target random non-zero bit string.
 
    Returns:
      str: A random bit string of length n that is not all zeroes.
    """
    # Generate random bits
    bits = [random.choice(['0', '1']) for _ in range(n)]

    # If all zeroes, set a random bit to 1
    if '1' not in bits:
      random_position = random.randint(0, n-1)
      bits[random_position] = '1'

    return ''.join(bits)

  @staticmethod
  def _xor_random_subset(bit_string: str, r: str):
    """
    Computes the XOR of a random subset of bits from the input bit string.

    Args:
      bit_string (str): A string containing only '0' and '1' characters.
      r (str): A bit string representing the random subset
        where each '1' represents the corresponding bit in bit_string is included in the xor

    Returns:
      str: The result of XORing the randomly selected bits ('0' or '1').

    Raises:
      ValueError: If the input string contains characters other than '0' or '1'.
    """
    # Validate input
    if not all(bit in '01' for bit in bit_string):
      raise ValueError("Input must be a bit string containing only '0' and '1' characters")
  
    if not bit_string:
      return '0'

    # Convert bit string to list of integers
    bits = [int(bit) for bit in bit_string]

    # Decide which bits to include
    selected_indices = [i for i, bit in enumerate(r) if bit == '1']

    # Compute XOR of selected bits
    result = 0
    for index in selected_indices:
      result ^= bits[index]  # XOR operation
  
    return str(result)

  @staticmethod
  def _concatenate_exp_with_bit(g: int, p: int, input_bit: str, x_binary: str):
    """
    Computes g^x mod p in binary format and concatenates it with the input bit.
    x is provided as a binary string.
    
    Args:
      g (int): Generator of the multiplicative group Z_p*.
      p (int): Prime modulus.
      input_bit (str): A single bit ('0' or '1') to concatenate.
      x_binary (str): Binary representation of the exponent x.
        
    Returns:
      - str: Concatenated bit string (binary(g^x mod p) + input_bit) 
        padded with zeroes to make the total length n + 1
        
    Raises:
      ValueError: If input_bit is not '0' or '1'.
      ValueError: If x_binary contains characters other than '0' or '1'.
    """
    # Validate input_bit
    if input_bit not in ['0', '1']:
      raise ValueError("Input bit must be either '0' or '1'")
 
    # Validate x_binary
    if not all(bit in '01' for bit in x_binary):
      raise ValueError("Input string must be a binary string containing only '0' and '1' characters")

    # Convert x from binary to integer
    x = int(x_binary, 2)

    # Compute g^x mod p
    result = pow(g, x, p)

    # Convert result to binary (removing '0b' prefix) and pad zeroes to make the total length n
    binary_result = bin(result)[2:].zfill(len(x_binary))

    # Concatenate with input_bit
    concatenated = binary_result + input_bit
    
    return concatenated

  def __init__(self, security_parameter: str, key: str):
    """
    Initializes the PRF by generating an n-bit safe prime and a generator and 
      initializing the security parameter and the key.

    Args:
      security_parameter (str): The security parameter input as a unary string.
      key (str): A key of same length as the security_parameter

    Raises:
      ValueError: 
        - If security_parameter is not a unary string,
        - If key is not a binary string,
        - If key and security_parameter are different length strings,
        - If security_parameter length < 3 or
        - If no safe prime is found after the maximum number of attempts in the initialization process
    """
    if not all(bit in '1' for bit in security_parameter):
      raise ValueError("Security Parameter must be a unary string")
    if not all(bit in '01' for bit in key):
      raise ValueError("Key must be a binary string")
    if len(key) != len(security_parameter):
      raise ValueError(f"Key length must match the input security parameter length {n}")
    n = len(security_parameter)
    if n < 3:
      raise ValueError("Security parameter must be at least 3 to generate a safe prime")

    self.n = n
    PRF._parameter_cache = PRF._load_cache()
    if n in PRF._parameter_cache:
      self._p, self._q, self._g = PRF._parameter_cache[n]
    else:
      self._p, self._q, self._g = self._generate_safe_prime_with_generator(n)
      PRF._parameter_cache[n] = (self._p, self._q, self._g)
      PRF._save_cache()
    self._key = key
    self._r = key
  
    if self._p is None:
      raise ValueError(f"Failed to Initialize - Please retry")

  def _is_probably_prime(self, num, k = 40):
    """
    Runs Miller-Rabin primality test.

    Args:
      num (int): The number to test for primality.
      k (int): The number of rounds of testing to perform.

    Returns:
      bool: True if the number is probably prime, False if it's definitely composite.
    """
    if num == 2 or num == 3:
      return True
    if num <= 1 or num % 2 == 0:
      return False

    # Write n - 1 as 2^r * d where d is odd
    r, d = 0, num - 1
    while d % 2 == 0:
      r += 1
      d //= 2

    # Witness loop
    for _ in range(k):
      a = random.randint(2, num - 2)
      x = pow(a, d, num)
      if x == 1 or x == num - 1:
        continue
      for _ in range(r - 1):
        x = pow(x, 2, num)
        if x == num - 1:
          break
      else:
        return False
    return True

  def _find_generator(self, p, q, t):
    """
    Finds a generator for Z_p* where p = 2q + 1 and q is prime by running at most t iterations
    Args:
      p (int): The prime number p = 2q + 1 whose generator is to be found.
      q (int): The Sophie Germain prime q where p = 2q + 1
      t (int): Maximum number of iterations allowed

    Returns:
      int or None: Returns a generator g of p if found or None if not found.
    """
    for _ in range(int(t)):
      g = random.randint(2, p - 2)

      # Check if g^q ≡ 1 (mod p)
      if pow(g, q, p) == 1:
        continue

      # Check if g^2 ≡ 1 (mod p)
      if pow(g, 2, p) == 1:
        continue
 
      # If g passes both tests, it's a generator
      return g
    return None
  
  def _generate_safe_prime_with_generator(self, n):
    """
    Generates an n-bit safe prime p, where p = 2q + 1 and q is also prime,
    along with a generator g for the multiplicative group Z_p*.
    
    Returns:
      tuple or (None, None, None): A tuple (p, q, g) or (None, None, None) if not found.
    """
    lower_bound = 2 ** (n - 1)
    upper_bound = 2 ** n - 1

    q_lower_bound = (lower_bound - 1) // 2
    q_upper_bound = (upper_bound - 1) // 2

    # Calculate maximum number of attempts: t = n³
    max_attempts = n ** 3

    attempts = 0
    while attempts < max_attempts:
      attempts += 1

      # Generate a random odd number for q within the range
      q_candidate = random.randrange(q_lower_bound, q_upper_bound + 1)
      if q_candidate % 2 == 0:
        q_candidate += 1
        if q_candidate > q_upper_bound:
          q_candidate = q_lower_bound + 1 if q_lower_bound % 2 == 0 else q_lower_bound

      # Check if q is prime
      if not self._is_probably_prime(q_candidate):
        continue

      # Calculate p = 2q + 1
      p_candidate = 2 * q_candidate + 1

      # Verify p is n-bit
      bit_length = p_candidate.bit_length()
      if bit_length != n:
        continue

      # Check if p is prime
      if self._is_probably_prime(p_candidate):
        # Found a safe prime p and its corresponding q
        # Polynomial iterations
        K = 2 ** n / (p_candidate - 1)
        # Now find a generator
        generator = self._find_generator(p_candidate, q_candidate, n * K)
        return (p_candidate, q_candidate, generator)

    # Return None if no safe prime is found after the maximum number of attempts
    return (None, None, None)

  def _get_params(self):
    """
    Return the generated parameters.
    
    Returns:
      tuple: A tuple (p, q, g) containing the safe prime, Sophie Germain prime, and generator.
    """
    return (self._p, self._q, self._g)

  def _verify_generator(self):
    """
    Verify that g is indeed a generator of the multiplicative group Z_p*.
    
    Returns:
      bool: True if g is a generator, False otherwise.
    """
    # Check that g^q ≢ 1 (mod p)
    if pow(self.g, self.q, self.p) == 1:
      return False

    # Check that g^2 ≢ 1 (mod p)
    if pow(self.g, 2, self.p) == 1:
      return False

    # Check that g^(p-1) ≡ 1 (mod p) (Fermat's Little Theorem)
    if pow(self.g, self.p - 1, self.p) != 1:
      return False
    return True

  def _base_prg(self, s: str):
    """
    This is the expanded PRG with expansion factor n + 1 that generates n + 1 pseudorandom bits 
      by taking an n bit seed as input
    Ref. Section 8.4.1 - Introduction to Modern Cryptography (Yehuda Lindell et al.) 3rd edition

    Args:
      s (str): Input n-bit string

    Returns:
      (str): A string of length n + 1 with pseudorandom bits

    Raises:
      ValueError: If the message is not a bit string or the input string length is not n
    """
    # Validate input
    if not all(bit in '01' for bit in s):
      raise ValueError("Input must be a bit string containing only '0' and '1' characters")
    if not len(s) == self.n:
      raise ValueError("Input must be a bit string of length n")

    if not s:
      return s

    # compute xor of random subset (determined by the key) of bits of s
    last_bit = self._xor_random_subset(s, self._r)
    assert(len(last_bit) == 1) # output should be 1 bit

    # compute g^x mod p and concatenate it with the xor from above step
    output = self._concatenate_exp_with_bit(self._g, self._p, last_bit, s)
    return output

  def _expanded_prg(self, s: str, p_n: int):
    """
    This is the expanded PRG with expansion factor n + p_n using the base PRG 
      that generates n + p_n pseudorandom bits by taking an n bit seed as input
    Ref. Section 8.4.2 - Introduction to Modern Cryptography (Yehuda Lindell et al.) 3rd edition

    Args:
      s (str): Input n-bit seed
      p_n (str): Expansion factor = n + p_n

    Returns:
      (str): A string of length n + p_n with pseudorandom bits

    Raises:
      ValueError: If the message is not a bit string or the input string length is not n
    """
    # Validate input
    if not all(bit in '01' for bit in s):
      raise ValueError("Input must be a bit string containing only '0' and '1' characters")
    if not len(s) == self.n:
      raise ValueError("Input must be a bit string of length n")

    if not s:
      return s

    output = s
    for _ in range(p_n):
      first_n = output[:self.n]
      remaining = output[self.n:]
      output = self._base_prg(first_n) + remaining

    return output

  def _pseudo_random_function(self, k: str, s: str):
    """
    This is the pseudorandom function based on the expanded PRG described above.
    Ref. Section 8.5 - Introduction to Modern Cryptography (Yehuda Lindell et al.) 3rd edition.

    Args:
      k (str): Input n-bit key
      s (str): Input n-bit string

    Returns:
      (str): A string of length n which is equal to F_k(s)

    Raises:
      ValueError: If the message is not a bit string or the input string length is not n
    """
    # Validate input
    if not all(bit in '01' for bit in s):
      raise ValueError("Input must be a bit string containing only '0' and '1' characters")
    if not len(s) == self.n:
      raise ValueError("Input must be a bit string of length n")

    output = k
    for i, bit in enumerate(s):
      if bit == '1':
        output = self._expanded_prg(output, self.n)[self.n:]
      else:
        output = self._expanded_prg(output, self.n)[:self.n]

    return output

  def call(self, m: str):
    """Wrapper exposing the pseudorandom function"""
    return self._pseudo_random_function(self._key, m)

  def __str__(self):
    """String representation of the parameters."""
    return f"PRF(n={self.n}, k={self._key} p={self._p}, q={self._q}, g={self._g})"

class MAC:
  """
  A class that takes as input a key and creates a MAC based on the pseudorandom function PRF
  """
  def __init__(self, security_parameter: str, key: str):
    """
    Initializes the MAC.

    Args:
      security_parameter (str): Security parameter as a unary n-bit string
      key (str): Input key with the same length as the security_parameter

    Raises:
      ValueError: If the message is not a bit string.
    """
    self._prf = PRF(security_parameter, key)
    self._key = key
    self.n = len(security_parameter)

  def generate_tag(self, message):
    """
    Generates a tag for the message.

    Args:
      message (str): The message for which the tag has to be generated.

    Returns:
      str: The generated tag.

    Raises:
      ValueError: If the message is not a bit string or length of the bit string is not 2(n - 1)
    """
    # Validate input
    if not all(bit in '01' for bit in message):
      raise ValueError("Message must be a bit string")
    if not len(message) == 2 * self.n - 2:
      raise ValueError("Invalid message length")

    m1 = message[:self.n - 1]
    m2 = message[self.n - 1:]

    # tag = F_k(0 || m1) || F_k(1 || m2)
    tag = self._prf.call('0' + m1) + self._prf.call('1' + m2)
    return tag

class Adversary:
  """
  A class that represents a PPT adversary that will attempt to attack the MAC.
  It takes as input a security parameter, a MAC oracle and a callback function 
  that will be used by the adversary to output a message and tag prediction
  """
  def __init__(self, security_parameter: str, oracle: Callable[[str], str], callback):
    """
    Initializes the relevant parameters

    Args:
      security_parameter (str): The security parameter for the MAC
      oracle (Callable[[str], str]): The oracle function that generates tags for messages
      callback (Callable[[str, str], None]): A callback function that takes the message and the tag 
        output by the adversary as input and verifies using canonical verification to check 
          - If the (message, tag) tuple was not previously queried to the oracle by the adversary and 
          - If the tag matches the tag output by the MAC for the same message
    """
    self.n = len(security_parameter)
    self.oracle = oracle
    self.callback = callback
  
  def attack(self):
    """
    Attacks the MAC.
    Strategy used for attack:
      1. Get tag t1 for 0 * (n - 2) || 1 * (n) using MAC oracle
      2. Get tag t2 for 1 * (n - 1)|| 0 * (n - 1) using MAC oracle
      3. Output a new message m = 0 * (n - 2) || 1 || 0 * (n - 1)
      4. Output tag t = t1[1...n] || t2[n + 1...2n]

      This attack is successful since the MAC is based on concatenation of output of two calls to 
        a deterministic pseudorandom function.
      The attacker can retrieve the full tag by masking half of the message for each query
        to get the tag for each half of the output.
    """
    first_message = '0' * (self.n - 2) + '1' + '1' * (self.n - 1)
    tag_1 = self.oracle(first_message)

    second_message = '0' * (self.n - 1) + '0' * (self.n - 1)
    tag_2 = self.oracle(second_message)

    output_message = '0' * (self.n - 2) + '1' + '0' * (self.n - 1)
    output_tag = tag_1[:n] + tag_2[n:]

    self.callback(output_message, output_tag)
    

class MAC_FORGE:
  """
  A class that takes as input a security parameter and a key, 
    creates an oracle based on the MAC and runs MAC-FORGE experiment
  """
  def __init__(self, security_parameter: str, key: str):
    """
    Initializes the relevant parameters.

    Args:
      security_parameter (str): Security parameter as a unary n-bit string
      key (str): Input key with the same length as the security_parameter
    """
    self.n = len(security_parameter)
    self._mac = MAC(security_parameter=security_parameter, key=key)
    self._Q = []
    self._result = None
    self._adversary_output = None

    def callback(m, tag):
      """
      This is the callback function for adversary to output a message and tag prediction
      It verifies 
        - the accuracy of the output tag
        - that the tag (m, tag) was not queried before by the adversary
      and sets the result of the experiment
      Args:
        m (str): The message
        tag (str): The tag
      """
      self._adversary_output = (m, tag)

      # Verify that the tuple of the message and the tag was never queried before
      for m_prime, tag_prime in self._Q:
        if m_prime != m or tag != tag_prime:
          continue
        self._result = 0
        return

      # Verify if the tag is accurate
      tag_prime = self._mac.generate_tag(m)
      if (tag == tag_prime):
        self._result = 1
      else:
        self._result = 0

    self._adversary = Adversary(security_parameter, oracle=self.oracle, callback=callback)

  def oracle(self, m: str, isAdversary: bool = True):
    """
    Generates a tag for the message and records the query by the adversary.

    Args:
      m (str): The message for which the tag has to be generated.

    Returns:
      str: The generated tag.
    """
    if not all(bit in '01' for bit in m):
      raise ValueError("Message must be a bit string")
    if not len(m) == 2 * self.n - 2:
      raise ValueError("Invalid message length")

    # generate tag
    tag = self._mac.generate_tag(m)
    # record the query by the adversary
    if isAdversary:
      self._Q.append((m, tag))
    return tag

  def run_experiment(self):
    """
    Runs the experiment and prints detailed information about the results of the experiment.

    Returns:
      int: The result of the experiment (1 for success, 0 for failure)
    """
    print("\n" + "="*50)
    print("EXPERIMENT RESULTS")
    print("="*50)

    # Run the attack
    self._adversary.attack()

    # Print oracle queries
    print("\nADVERSARY QUERIES:")
    print("-"* (self.n * 5))
    if len(self._Q) == 0:
      print("No queries were made to the oracle.")
    else:
      print(f"{'#':<5} {'Message':<{2*self.n}} {'Tag':<{self.n}}")
      print("-"* (self.n * 5))
      for i, (message, tag) in enumerate(self._Q):
        print(f"{i+1:<5} {message:<{2*self.n}} {tag:<{self.n}}")

    # Print adversary output
    print("\nADVERSARY OUTPUT:")
    print("-"*(self.n * 3))
    if self._adversary_output is None:
      print("Adversary did not produce any output.")
    else:
      m, tag = self._adversary_output
      true_tag = self._mac.generate_tag(m)

      print(f"Message:     {m}")
      print(f"Forged tag:  {tag}")
      print(f"Correct tag: {true_tag}")
      print(f"Tag match:   {'Yes' if tag == true_tag else 'No'}")

      # Check if this was a query
      was_queried = any(m == m_q and tag == t_q for m_q, t_q in self._Q)
      print(f"Was queried: {'Yes' if was_queried else 'No'}")

    # Print experiment result
    print("\nEXPERIMENT RESULT:")
    print("-"*50)
    result_str = "SUCCESS" if self._result == 1 else "FAILURE"
    print(f"Result: {result_str} ({self._result})")

    # Print explanation
    print("\nEXPLANATION:")
    print("-"*50)
    if self._result == 1:
      print("The adversary successfully forged a valid tag for a message")
      print("that was not previously queried. The MAC is not secure.")
    else:
      if self._adversary_output is None:
        print("The adversary failed to produce any output.")
      elif any(m == self._adversary_output[0] and tag == self._adversary_output[1] 
              for m, tag in self._Q):
        print("The adversary returned a message-tag pair that was")
        print("previously queried. This does not constitute a forgery.")
      else:
        print("The adversary produced an invalid tag for the message.")
        print("The forgery attempt failed.")

    print("="*50 + "\n")

    return self._result

def generate_random_bit_string(n):
  """
  Generates a random binary string of length n.

  Args:
    n: Length of the binary string to generate

  Returns:
    A string of '0's and '1's of length n
  """
  return ''.join(random.choice(['0', '1']) for _ in range(n))

if __name__ == "__main__":
  n = None
  key = None
  mac = None

  print("\nMAC-FORGE Experiment")
  print("="*20)
  print()

  while True:
    # Only ask for configuration if not already set
    if mac is None:
      # Get parameter n
      n_choice = input("Enter security parameter n or press Enter for default (n = 30): ").strip()
      if n_choice:
        try:
          n = int(n_choice)
          if n <= 0:
            print("Error: n must be a positive integer.")
            continue
        except ValueError:
          print("Error: Please enter a valid integer for n.")
          continue
      else:
        n = 30

      # Get key
      key_choice = input(f"Enter a key of length {n} or press Enter for a random key: ").strip()
      if key_choice:
        if len(key_choice) != n:
          print(f"Error: Key must be exactly {n} characters long.")
          continue
        key = key_choice
      else:
        key = generate_random_bit_string(n)
        print(f"Generated random key: {key}")
      
      try:
        mac = MAC_FORGE(security_parameter='1' * n, key=key)
        print(f"MAC initialized with n = {n} and key of length {n}")
      except ValueError as e:
        print(f"Error initializing MAC: {e}")
        mac = None
        continue
    
    # Display menu
    print("\nOptions:")
    print("1. Run MAC-FORGE experiment")
    print("2. Interact with the MAC Oracle")
    print("3. Reset configuration")
    print("4. Quit the program")
    
    choice = input("\nEnter your choice (1-4): ").strip()
  
    if choice == "1":
      try:
        print("\nRunning MAC-FORGE experiment...")
        mac.run_experiment()
      except Exception as e:
        print(f"Error during experiment: {e}")
    
    elif choice == "2":
      print("\nMAC Oracle Interaction")
      print("=====================")
      message = input("Enter a message to generate MAC (or 'back' to return to menu): ")
      
      while message.lower() != 'back':
        try:
          tag = mac.oracle(message, False)
          print(f"Message: {message}")
          print(f"MAC tag: {tag}")
        except Exception as e:
          print(f"Error generating MAC: {e}")
        
        message = input("\nEnter another message (or 'back' to return to menu): ")
    
    elif choice == "3":
      print("\nResetting configuration...")
      n = None
      key = None
      mac = None

    elif choice == "4":
      print("\nExiting MAC-FORGE program. Goodbye!")
      break

    else:
      print("\nInvalid choice. Please enter a number between 1 and 4.")
