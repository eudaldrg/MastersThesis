#define BOOST_TEST_MODULE UtilsUnitTests

#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_parameters.hpp>

using namespace boost::unit_test;

struct Init {
    Init()
    {
        boost::unit_test::unit_test_log.set_threshold_level( boost::unit_test::log_messages);
//        if (framework::master_test_suite().argc == 2 && std::string(framework::master_test_suite().argv[1]) == "--suppress_xan_logs")
//            LogNS::SetDefaultLogL(LogL::UnusedLevel);
//        else if(  runtime_config::get<log_level>(runtime_config::btrt_log_level) <= log_messages )
//            LogNS::SetDefaultLogL(LogL::Debug4);
//        else
//            LogNS::SetDefaultLogL(LogL::Error);
    }
    ~Init()  = default;
};

BOOST_GLOBAL_FIXTURE( Init );
