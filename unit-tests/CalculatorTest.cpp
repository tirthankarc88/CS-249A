#include "gtest/gtest.h"
#include "Calculator.h"

TEST(Calculator, addsTwoNewNumbers) {
  Calculator calc;
  calc.firstOperandIs(2);
  calc.secondOperandIs(5);
  ASSERT_TRUE(7 == calc.add());
}

TEST(Calculator, multipliesTwoNewNumbers) {
  Calculator calc;
  calc.firstOperandIs(2);
  calc.secondOperandIs(5);
  ASSERT_TRUE(10 == calc.multiply());
}

TEST(Calculator, switchingOperandsForAdd) {
  Calculator calc;
  calc.firstOperandIs(2);
  calc.secondOperandIs(5);
  calc.secondOperandIs(9);
  ASSERT_TRUE(11 == calc.add());
}

TEST(Calculator, switchingOperandsForMultiply) {
  Calculator calc;
  calc.firstOperandIs(2);
  calc.secondOperandIs(5);
  calc.secondOperandIs(9);
  ASSERT_TRUE(18 == calc.multiply());
}
