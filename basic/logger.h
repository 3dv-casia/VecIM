#ifndef CM_LOGGER_H
#define CM_LOGGER_H

#include <string>
#include <fstream>

#include <boost/log/trivial.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/sources/global_logger_storage.hpp>
#include <boost/log/sinks/sync_frontend.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>

namespace cm {

typedef boost::log::trivial::severity_level severity_level;

/*!
 * \brief Global narrow-char thread-safe logger.
 */
BOOST_LOG_INLINE_GLOBAL_LOGGER_DEFAULT(Logger, boost::log::sources::severity_logger_mt<severity_level>)

#define LOG_TRACE BOOST_LOG_SEV(cm::Logger::get(), cm::severity_level::trace)
#define LOG_DEBUG BOOST_LOG_SEV(cm::Logger::get(), cm::severity_level::debug)
#define LOG_INFO BOOST_LOG_SEV(cm::Logger::get(), cm::severity_level::info)
#define LOG_WARNING BOOST_LOG_SEV(cm::Logger::get(), cm::severity_level::warning)
#define LOG_ERROR BOOST_LOG_SEV(cm::Logger::get(), cm::severity_level::error)
#define LOG_FATAL BOOST_LOG_SEV(cm::Logger::get(), cm::severity_level::fatal)

/*!
 * \brief Custom colored formater.
 * \param rec boost log record_view
 * \param strm boost log formatting ostream
 */
inline void coloring_formatter(
  boost::log::record_view const &rec,
  boost::log::formatting_ostream &strm) {
  // color w.r.t severity level
  auto severity = rec[boost::log::trivial::severity];
  if (severity) {
    // set the leading color
    switch (severity.get()) {
    case severity_level::trace:
      strm << "[TRACE]: ";
      break;
    case severity_level::debug:
      strm << "\033[34m[DEBUG]\033[0m: ";
      break;
    case severity_level::info:
      strm << "\033[32m[INFO]\033[0m: ";
      break;
    case severity_level::warning:
      strm << "\033[33m[WARNING]\033[0m: ";
      break;
    case severity_level::error:
      strm << "\033[31m[ERROR]\033[0m: ";
      break;
    case severity_level::fatal:
      strm << "\033[41m[FATAL]\033[0m: ";
      break;
    default:
      break;
    }
  }

  // format the message here
  strm << rec[boost::log::expressions::smessage];
}

/*!
 * \brief Initializes boost logger directory and format.
 * \param level severity level
 * \param log_path log file path, no file log if it is empty
 * \param ansi_color set true to turn on console color
 */
inline void initialize_logger(
  const severity_level &level = severity_level::info,
  const std::string &log_path = "",
  const bool ansi_color = true) {
  // namespace
  namespace log = boost::log;
  namespace sinks = boost::log::sinks;
  namespace keywords = boost::log::keywords;

  log::register_simple_formatter_factory<severity_level, char>("Severity");
  if (!log_path.empty()) {
    log::add_file_log(
      keywords::file_name = log_path,
      keywords::rotation_size = 10 * 1024 * 1024,
      keywords::time_based_rotation = sinks::file::rotation_at_time_point(0, 0, 0),
      keywords::format = "[%TimeStamp%][%ThreadID%][%Severity%]: %Message%"
    );
  }
  if (ansi_color) {
    // construct the sink
    typedef sinks::synchronous_sink<sinks::text_ostream_backend> text_sink;
    boost::shared_ptr<text_sink> sink = boost::make_shared<text_sink>();

    // text_ostream_backend supports adding several streams.
    boost::shared_ptr<std::ostream> stream(&std::clog, boost::null_deleter());
    sink->locked_backend()->add_stream(stream);

    // sink color formatter
    sink->set_formatter(&coloring_formatter);

    // register the sink in the log core
    log::core::get()->add_sink(sink);
  }
  else
    log::add_console_log(std::cout, keywords::format = "[%Severity%]: %Message%");

  log::add_common_attributes();
  log::core::get()->set_filter(log::trivial::severity >= level);
}

/*!
 * \brief Change serevity level during runtime.
 * \param level severity level
 */
inline void set_severity_level(const severity_level &level) {
  boost::log::core::get()->set_filter(boost::log::trivial::severity >= level);
}

} // namespace cm

#endif // CM_LOGGER_H
