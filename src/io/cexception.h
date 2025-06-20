#ifndef CEXCEPTION_H
#define CEXCEPTION_H

#include <string>
#include <iostream>
#include <exception>

class CException : public std::exception
{
    public:
        CException(std::string Message):message(Message) {}
        virtual ~CException() throw() {}
        std::string report() {return message;}
        virtual const char* what() const throw() {
            std::string tmp("\n\nERROR: ");
            tmp+=message+std::string("\n\n");
            return tmp.c_str();
        }
    protected:
        std::string message;
    private:
};

#endif // CEXCEPTION_H
