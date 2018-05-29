#ifndef _CMD_EXCEPTION_H_
#define _CMD_EXCEPTION_H_

#include <exception>
#include <string>

#define CMD_EXCEPTION(msg) CommandException( std::string(msg) )

class CommandException : public std::exception {
public:
	explicit CommandException(const std::string& _message) {
		this->str_msg_error = _message;
	}

	char const* what() const override {
		return str_msg_error.c_str();
	}

private:
	std::string str_msg_error;
};

#endif // _CMD_EXCEPTION_H_
