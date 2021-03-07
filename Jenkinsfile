pipeline {
  parameters{string(defaultValue
                    : '', description
                    : 'api-token', name
                    : 'API_TOKEN') string(defaultValue
                                          : '', description
                                          : 'Token for readthedocs', name
                                          : 'RTD_TOKEN')
                 string(defaultValue
                        : '', description
                        : 'buildable phid', name
                        : 'BUILD_TARGET_PHID') string(defaultValue
                                                      : '', description
                                                      : 'Commit id', name
                                                      : 'COMMIT_ID')
                     string(defaultValue
                            : '', description
                            : 'Diff id', name
                            : 'DIFF_ID')
                         string(defaultValue
                                : 'PHID-PROJ-qr4btg2lycr6fpvkkmb5', description
                                : 'ID of the project', name
                                : 'PROJECT_ID')}

  environment{
    PHABRICATOR_HOST = 'https://c4science.ch/api/' PYTHONPATH = sh
    returnStdout: true,
    script: 'echo ${WORKSPACE}/tests/ci/script/'
  }

  agent {
    dockerfile { additionalBuildArgs '--tag cracklet-environment' }
  }

  stages {
    stage('SCM Checkout') {
      steps {
        checkout scm: [$class:'GitSCM',
                                    branches:scm.branches,
                                  extensions:[[$class:'SubmoduleOption',
                                                 recursiveSubmodules:true, ]],
                           userRemoteConfigs:scm.userRemoteConfigs]
      }
    }

    stage('Configure') {
      steps{
        sh ""
           "#!/bin/bash
            set -
            o pipefail mkdir - p build cd build cmake -
            DCRACKLET_PYTHON_INTERFACE : BOOL =
            TRUE - DCRACKLET_EXAMPLES : BOOL =
                TRUE - DCRACKLET_TESTS : BOOL =
                    TRUE - DCMAKE_BUILD_TYPE : STRING =
                        RelWithDebInfo..| tee../ configure.txt ""
                                                               "
      } post {
        failure {
          uploadArtifact('build/configure.txt', 'Configure') sh ""
                                                                "
              rm -
              rf build ""
                       "
        }
      }
    }

    stage('compile') {
      steps{sh '''#!/bin/bash set - o pipefail make - C build / src |
            tee build / compilation.txt
		      ''' } post {
        failure { uploadArtifact('build/compilation.txt', 'Compilation') }
      }
    }

    stage('compile python') {
      steps{sh '''#!/bin/bash set -
                o pipefail

                    make -
                C build / python |
            tee build / compilation_python.txt
			     ''' } post {
        failure {
          uploadArtifact('build/compilation_python.txt', 'Compilation_Python')
        }
      }
    }
  }
}