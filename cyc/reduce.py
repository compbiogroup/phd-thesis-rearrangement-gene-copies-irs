import subprocess

'''
Automatic reduction of instances for testing
'''


with open("in.txt", "r") as file:
    s = file.readline().split()
    # file.readline()
    p = file.readline().split()
    # file.readline()

for j in [1,2,3]:
    i = 1
    while i < (len(s) if j < 3 else len(p)) - 1:
        if (j == 1):
            v = s[i]
        if j > 1 or v in p:
            s_ = s.copy()
            p_ = p.copy()
            if j == 1:
                s_.pop(i)
                p_.remove(v)
            if j == 2:
                s_.pop(i)
            if j == 3:
                p_.pop(i)

            with open("in_tmp.txt", "w") as file:
                print(" ".join(s_), file = file)
                # print(" ".join(["0" for _ in range(len(s_) - 1)]), file = file)
                print(" ".join(p_), file = file)
                # print(" ".join(["0" for _ in range(len(p_) - 1)]), file = file)

            # timeoutSeconds = 45
            # try:
            #     subprocess.run(['./MCISP2K', '-i', "in_tmp.txt", "-o", "out.txt"], timeout = timeoutSeconds)
            # except subprocess.TimeoutExpired:
            #     print("timeout")
            try:
                # subprocess.run(['./dec', "ga", '-i', "in_tmp.txt", "-k", "100", "-o", "out"], check=True)
                subprocess.run(['./dist', "reversal_noir", '-i', "in_tmp.txt", "-k", "1", "-o", "out", "-z", "-e"], check=True)
            except subprocess.CalledProcessError:
                print("error")
                s = s_
                p = p_
                print(s)
                print(p)
                i -= 1
            else:
                print("solved")
        i += 1

    with open("in_" + str(j) + ".txt", "w") as file:
        print(" ".join(s), file = file)
        # print(" ".join(["0" for _ in range(len(s) - 1)]), file = file)
        print(" ".join(p), file = file)
        # print(" ".join(["0" for _ in range(len(p) - 1)]), file = file)
