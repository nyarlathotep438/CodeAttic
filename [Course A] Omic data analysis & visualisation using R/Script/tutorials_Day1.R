#Task3
X = 5
Y = 10
XY = X + Y

#Task4
X = 1
XY = X + Y

#Task6
numbers1 = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)
numbers2 = c(22, 24, 26, 28, 30)
numbers = c(numbers1, numbers2)

#Task7
numbers_mean = mean(numbers)

#Task8
numbers[3]
index12 = numbers[12]
numbers[666]

numbers[5] = 999
numbers[6] = "potato"
numbers[7] = TRUE

#Task9
vector_for_loop = c(101:200)
for (index in 1:100)
{
  print(index)
}

for (index in 1:100)
{
  print(vector_for_loop[index])
}

#Task10
x = 10

x > 0

x < 0

x != 0

x == 0

if (x > 0) {print("I am bigger than zero")}


if (x < 0) {print("I am smaller than zero")}


if (x != 0) {print("I am not zero")}


if (x == 0) {print("I am zero")}

#Task11
for (number in 1:100)
{
  if (number > 50){
    print(number)
  }
}

9 %% 2 # will return 1
8 %% 2 # will return 0
7 %% 2 # will return 1
6 %% 2 # will return 0