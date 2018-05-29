#ifndef __COMMAND_MANAGER_H__
#define __COMMAND_MANAGER_H__

#include <stack>
#include <vector>
#include "Exception.h"

class CommandManager; // Predestination

class Command {
public:
	explicit Command(std::vector<int32_t>& _vector) : vector(_vector) {}
	virtual void undo() = 0;
	virtual void redo() = 0;
	virtual void execute() = 0;
protected:
	virtual void doExecute() = 0;
	std::vector<int32_t>& vector;

	friend CommandManager;
};

class CommandManager {
public:
	static CommandManager* getInstance() {
		static CommandManager* instance;
		if (instance == nullptr)
			instance = new CommandManager();
		return instance;
	}

	void cmd_execute(Command* cmd) {
		cmdStackUndo.push(cmd);
		cmd->doExecute();
	}

	void cmd_undo() {
		if (cmdStackUndo.size() <= 0)
			throw CMD_EXCEPTION("Error in command Undo. Stack is Empty");

		cmdStackUndo.top()->undo();
		cmdStackRedo.push(cmdStackUndo.top());
		cmdStackUndo.pop();
	}

	void cmd_redo() {
		if(cmdStackUndo.size() <= 0)
			throw CMD_EXCEPTION("Error in command Redu. Stack is Empty");

		cmdStackRedo.top()->redo();
		cmdStackUndo.push(cmdStackRedo.top());
		cmdStackRedo.pop();
	}

private:
	CommandManager() {}
	~CommandManager() {}
	CommandManager(CommandManager const&) = delete;
	CommandManager& operator= (CommandManager const&) = delete;

	std::stack<Command*> cmdStackUndo;
	std::stack<Command*> cmdStackRedo;
};

class CMDInsertBackElement : public Command {
public:
	CMDInsertBackElement(std::vector<int32_t>& vector, int value) : Command(vector) { this->value = value; }
	void execute() override { CommandManager::getInstance()->cmd_execute(this); }
private:
	int32_t value;
	void doExecute() override { vector.push_back(value); }
	void redo() override { doExecute(); }
	void undo() override { if (!vector.empty()) vector.pop_back(); }
};

class CMDDeleteBackElement : public Command {
public:
	CMDDeleteBackElement(std::vector<int>& vector) : Command(vector) {
		if (vector.empty())
			throw CMD_EXCEPTION("Error in command DeleteBackElement. Vector is Empty");

			saveValue = vector.back();
	}

	void execute() override { CommandManager::getInstance()->cmd_execute(this); }

private:
	int32_t saveValue;
	void doExecute() override {
		if (vector.empty())
			throw CMD_EXCEPTION("Error execute command DeleteBackElement. Vector is Empty");

		vector.pop_back();
	}
	void undo() override { vector.push_back(saveValue); }
	void redo() override { doExecute(); }
};

class CMDSetElement : public Command {
public:
	explicit CMDSetElement(std::vector<int32_t>& vector, int32_t value, int32_t index) : Command(vector) {
		if (vector.size() <= static_cast<size_t>(index))
			throw CMD_EXCEPTION("Error in command SetElement. Invalid index");

		this->value = value;
		this->saveIndex = index;
		this->saveValue = vector[index];
	}
	void execute() override { CommandManager::getInstance()->cmd_execute(this); }

private:
	int32_t saveValue, saveIndex, value;

	void doExecute() override {
		if (vector.size() <= static_cast<size_t>(saveIndex))
			throw CMD_EXCEPTION("Error execute command SetElement. Invalid index");

		vector[saveIndex] = value;
	}

	void undo() override {
		if (vector.size() > static_cast<size_t>(saveIndex))
			vector[saveIndex] = saveValue;
	}

	void redo() override { doExecute(); }
};

#endif // __COMMAND_MANAGER_H__
