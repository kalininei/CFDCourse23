@page impl-details Детали реализации
@tableofcontents

# Основные классы {#MainClasses}

tal::Executor
--- *src/tal_executor/tal_executor.hpp*

Движок, через который осуществляется взаимодействие с пользователем

Содержит в себе классы:

tal::ModulesContainer
--- *src/tal_executor/detail/tal_modules_container.hpp*

Контейнер, хранящий вычислительные модули. Поскольку на одну вычислительную роль может быть определён только один вычислительный модуль, то рабочие и входные модули хранятся в `std::tuple<>` - структуре.
Аггрегаторы, которых может быть произвольное количество, хранятся в векторе.

tal::FlowGraph
--- *src/tal_executor/detail/tal_flow_graph.hpp*

Дерево исполнения, в котором реализованы алгоритмы связывания и запуска графа tbb

Исполнение графа производится на платформе `oneapi/tbb/flow`.
Сборка графа происходит в классе `tal::FlowGraph`, однако основной
функционал для связывания и исполнения дерева зашит в модули.

# Схема дерева TBB для исполнения рабочего модуля {#TbbNodes}

Граф TBB, получающийся в результате вызова `tal::Executor::assemble()`, имеет следующую структуру. В красном прямоугольнике - узлы, отвечающие за вызов одного рабочего модуля.
@image html tbb_nodes.png width=400

Входные данные генерируются узлом `input_node`.

Рабочая функция запускается из узла `function_node`. Запускается в последовательном режиме `tbb::flow::serial`.
То есть не допускается одновременного исполнения одной рабочей функции для разных входных данных.
Тем не менее рабочие функции разных модулей могут исполняться асинхронно, если это позволяет схема связывания.

Узел `join_node` собирает аргументы в кортеж.
Имеет tbb::flow::reserving политику и требует внешний буфферный узел для накопления данных.

`buffer` в зависимости от контекста 
[TOC]
- либо `buffer_node` для обычного связывания
- либо `write_once_node` для опционального (обёртка `tal::Optional`) связывания в случае, когда данные на этот вход отсутствуют.
  Тогда этот узел будет всегда возвращать `nullptr`.

# Реализация класса модулей и ролей {#Modules}
Диаграмма наследования класса IModule. Синие прямоугольники - абстрактные классы, определяемые в библиотеке.
Красные -- классы, реализуемые пользователем.
@image html modules.png width=800

**tal::IModule** 
--- *src/tal_executor/detail/tal_imodule.hpp*

Базовый нешаблонный класс предоставляет
- общий для всех модулей интерфейс,
- логгер, доступный по вызову logger()
- механизм регистрации модулей по текстовому имени, необходимый для сериализации/десериализации модулей

**tal::SenderModule**
--- *src/tal_executor/detail/tal_sender_module.hpp*

На этом уровне известен тип возвращаемого значения, определённый шаблонным параметром `Role`.

Предоставляет
- псевдонимы для возвращаемого значения и указателя для него
- алгоритмы для работы с дампами
- доступ к последнему успешному результату

**tal::ReceiverModule**
--- *src/tal_executor/detail/tal_receiver_module.hpp*

- хранит информацию по одной входной зависимости (зависимость от нескольких аргументов реализована через наследование от нескольких `ReceiverModule`)
- предоставляет псевдонимы для принимаемого значения и указателя для него
- хранит буферный узел tbb (`buffer_node` или `write_once_node`) и реализует внешнее связывание модуля в дереве tbb
- производит проверку входных данных

**tal::InvokerModule**
--- *src/tal_executor/detail/tal_invoker_module.hpp*
- обеспечивает запуск рабочей функции, а также запускает пре/пост процессоры вокруг неё в функции `run()`
- ведёт счётчик исполнения

**tal::DependentModule**
--- *src/tal_executor/detail/tal_dependent_module.hpp*
- обеспечивает связь со всеми `ReceiverModule`.
- собирает и хранит узел `tbb::flow::join_node`, применяемый для упаковки входных данных перед подачей в тело исполнения

**tal::IInputModule**
--- *src/tal_executor/detail/tal_iinput_module.hpp*
- Интерфейс для входных модулей (функция `activate()`)

**tal::BasicInputModule**
--- *src/tal_executor/tal_basic_input_module.hpp*
- Обеспечивает построение узла tbb::flow::input_node, вызывающего пользовательскую функцию `process()`

**tal::BasicWorkerModule**
--- *src/tal_executor/tal_basic_worker_module.hpp*
- Обеспечивает построение узла tbb::flow::function_node, вызывающего пользовательскую функцию `process()`

**tal::BasicAggregatorModule**
--- *src/tal_executor/tal_basic_aggregator_module.hpp*
- Обеспечивает построение узла tbb::flow::function_node, вызывающего пользовательскую функцию `aggregate()`