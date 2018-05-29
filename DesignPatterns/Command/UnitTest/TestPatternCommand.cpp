#include "stdafx.h"
#include "CppUnitTest.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTestPatternCommand
{		
	TEST_CLASS(VectorCommandTest)
	{
	public:
		VectorCommandTest() {
			size = 3;
			vector.clear();
			vector_expected.clear();
			vector.resize(size);
			vector_expected.resize(size);
			std::srand(unsigned(std::time(0)));
			std::generate(vector_expected.begin(), vector_expected.end(), [&]() { return std::rand(); });
			vector = vector_expected;
		}

		~VectorCommandTest() { }
		
		TEST_METHOD(test_deleteBackElement) {
			try {
				CMDDeleteBackElement* cmd = new CMDDeleteBackElement(vector);
				cmd->execute();

				Assert::IsTrue(vector != vector_expected, L"Error! Vectors are equal");
				Assert::IsTrue(vector.size() == vector_expected.size() - 1, L"Error! Does not correspond to the sizes");
			}
			catch (CommandException &e) {
				Logger::WriteMessage(e.what());
				Assert::Fail(L"Error! See the output");
			}
		}

		TEST_METHOD(test_insertBackElement) {
			value = 500;
			try {
				CMDInsertBackElement* cmd = new CMDInsertBackElement(vector, value);
				cmd->execute();

				Assert::IsTrue(vector.size() == vector_expected.size() + 1, L"Error! Does not correspond to the sizes");
				Assert::AreEqual(value, vector.back(), L"Error! Invalid value");
			}
			catch (CommandException &e) {
				Logger::WriteMessage(e.what());
				Assert::Fail(L"Error! See the output");
			}
		}

		TEST_METHOD(test_setElement) {
			value = 750, index = 1;
			try {
				CMDSetElement* cmd = new CMDSetElement(vector, value, index);
				cmd->execute();

				Assert::IsTrue(vector.size() == vector_expected.size(), L"Error! Does not correspond to the sizes");
				Assert::AreNotEqual(vector[index], vector_expected[index], L"Error! Invalid value. Is Equal");
			}
			catch (CommandException &e) {
				Logger::WriteMessage(e.what());
				Assert::Fail(L"Error! See the output");
			}
		}

		TEST_METHOD(test_undo) {
			value = 350, index = 1;
			try {
				CMDSetElement* cmd = new CMDSetElement(vector, value, index);
				cmd->execute();

				Assert::IsTrue(vector.size() == vector_expected.size(), L"Error! Does not correspond to the sizes");
				Assert::AreNotEqual(vector[index], vector_expected[index], L"Error! Invalid value. Is Equal");

				CommandManager::getInstance()->cmd_undo();

				Assert::AreEqual(vector[index], vector_expected[index], L"Error! Invalid value. Is Not Equal");
			}
			catch (CommandException &e) {
				Logger::WriteMessage(e.what());
				Assert::Fail(L"Error! See the output");
			}
		}

		TEST_METHOD(test_redo) {
			value = 350;
			try {
				CMDInsertBackElement* cmd = new CMDInsertBackElement(vector, value);
				cmd->execute();

				Assert::IsTrue(vector.size() == vector_expected.size() + 1, L"Error! Does not correspond to the sizes");
				Assert::AreEqual(value, vector.back(), L"Error! Invalid value");

				CommandManager::getInstance()->cmd_undo();

				Assert::IsTrue(vector.size() == vector_expected.size(), L"Error! Does not correspond to the sizes");
				Assert::AreEqual(vector.back(), vector_expected.back(),  L"Error! Invalid value");

				CommandManager::getInstance()->cmd_redo();

				Assert::IsTrue(vector.size() == vector_expected.size() + 1, L"Error! Does not correspond to the sizes");
				Assert::AreEqual(value, vector.back(), L"Error! Invalid value");
			}
			catch (CommandException &e) {
				Logger::WriteMessage(e.what());
				Assert::Fail(L"Error! See the output");
			}
		}

	private:
		int32_t size, index, value;
		std::vector<int32_t> vector;
		std::vector<int32_t> vector_expected;
	};
}