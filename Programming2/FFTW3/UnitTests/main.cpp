#define BOOST_TEST_MODULE UtilsUnitTests

#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_parameters.hpp>

using namespace boost::unit_test;

struct Init {
    Init()
    {
//        runtime_config::get<log_level>(runtime_config::btrt_log_level) = log_messages;
    }
    ~Init()  = default;
};

BOOST_GLOBAL_FIXTURE( Init );
